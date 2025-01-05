use super::{poseidon_parameters, solve_linear_system, FieldMT, FieldPath, Fp, F};
use ark_crypto_primitives::{crh::poseidon, CRHScheme};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, Polynomial};

// This function executes one folding step in the Generalized FRI algorithm
fn fold(f: &DensePolynomial<F>, l: u8, theta: F) -> DensePolynomial<F> {
    let mut g_polys: Vec<Vec<F>> = Vec::new();
    let d = (((Polynomial::degree(f)) / (l as usize)) as f32).floor() as usize + 1;
    for j in 0..l {
        let th = theta.pow(&[j as u64]);
        // The g are the g_i such that f(x) = x^i*g_i(x^l)
        let mut g: Vec<F> = vec![F::from(0); d];

        for (i, coeff) in f.coeffs.iter().enumerate() {
            if (j as u64 > i as u64) && ((j as u64 - i as u64) % l as u64 == 0) {
                g[i] += coeff * (&th);
            } else if (i as u64 >= j as u64) && ((i as u64 - j as u64) % l as u64 == 0) {
                g[(i - j as usize) / (l as usize)] += coeff * (&th);
            }
        }
        g_polys.push(g);
    }
    let mut final_g = vec![F::from(0); d];
    for j in 0..d {
        for poly in g_polys.iter() {
            final_g[j] += poly[j];
        }
    }
    while final_g.last().unwrap() == &F::from(0) {
        final_g.remove(final_g.len() - 1);
    }
    DensePolynomial { coeffs: final_g }
}

// Computes the Merkle tree of the folded f on the evaluation domain r<s>
fn round_commit(f_folded: &DensePolynomial<F>, s: &F, r: &F, s_ord: &u64) -> (FieldMT, Fp, Vec<F>) {
    let leaf_crh_params = poseidon_parameters();
    let two_to_one_params = leaf_crh_params.clone();

    let D: Vec<F> = (0..*s_ord).into_iter().map(|i| r * &s.pow([i])).collect();
    let point_vec: Vec<F> = D.into_iter().map(|x| f_folded.evaluate(&x)).collect();
    let k = (((*s_ord as f32).log2()).ceil()) as u32;
    let mut eval_vec: Vec<Vec<Fp>> = point_vec.iter().map(|x| vec![x.c0(), x.c1()]).collect();
    eval_vec = vec![eval_vec, vec![vec![Fp::from(0), Fp::from(0)]; 2u64.pow(k) as usize - *s_ord as usize]].concat();
    let eval_vec_slice: Vec<&[Fp]> = eval_vec.iter().map(|x| x.as_slice()).collect();

    let mtree: FieldMT = FieldMT::new(&leaf_crh_params, &two_to_one_params, eval_vec_slice).unwrap();
    (mtree.clone(), mtree.root(), point_vec)
}

// Returns (Merkle Trees, Merkle Roots, Evaluations)
fn commit(
    f: DensePolynomial<F>, l_list: Vec<usize>, mut s: F, mut r: F, mut s_ord: u64,
) -> (Vec<FieldMT>, Vec<Fp>, Vec<Vec<F>>) {
    let mut mtrees: Vec<FieldMT> = Vec::new();
    let mut points: Vec<Vec<F>> = Vec::new();
    let mut roots: Vec<Fp> = Vec::new();

    let (first_mt, first_root, first_points) = round_commit(&f, &s, &r, &s_ord);
    mtrees.push(first_mt);
    points.push(first_points);
    roots.push(first_root);
    let params = poseidon_parameters();
    let mut theta_vec: Vec<F> = Vec::new();
    theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0)));
    let n_rounds = l_list.len();
    let folded_polys: &mut Vec<DensePolynomial<F>> = &mut Vec::new();
    folded_polys.push(f);
    for i in 0..n_rounds {
        folded_polys.push(fold(folded_polys.last().unwrap(), l_list[i] as u8, theta_vec[i]));
        r = r.pow([l_list[i] as u64]);
        s = s.pow([l_list[i] as u64]);
        s_ord /= l_list[i] as u64;
        let (m, r, p) = round_commit(folded_polys.last().unwrap(), &s, &r, &s_ord);
        roots.push(r);
        mtrees.push(m);
        points.push(p);
        theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0)));
    }
    (mtrees, roots, points)
}

fn query_at_index(
    mt1: FieldMT, mt2: FieldMT, points1: Vec<F>, points2: Vec<F>, index: usize, l: usize, n: usize,
) -> (Vec<FieldPath>, Vec<F>) {
    let mut paths: Vec<FieldPath> = Vec::new();
    // Query MT2 at position index, save point to var, path to vec
    let path_main: FieldPath = mt2.generate_proof(index % (n / l)).unwrap();
    paths.push(path_main);
    let y: F = points2[index % (n / l)];
    // Query MT1 at positions (index +kn/l) % n, save points to vec, paths to vec paths
    let mut points: Vec<F> = vec![y];
    for k in 0..l {
        let i: usize = (index + (k * n) / l) % n;
        let path_current: FieldPath = mt1.generate_proof(i).unwrap();
        paths.push(path_current);
        points.push(points1[i]);
    }

    // Output paths to be verified later
    (paths, points)
}

fn query(
    mtrees: Vec<FieldMT>, points: Vec<Vec<F>>, l_list: Vec<usize>, n: usize, alpha: usize, grinding_param: u8,
) -> (Vec<FieldPath>, Vec<F>, Vec<usize>) {
    let mut paths: Vec<FieldPath> = Vec::new();
    let mut queried_points: Vec<F> = Vec::new();

    // Vec of indices to query consistency for PolyIOP
    let mut indices_first: Vec<usize> = Vec::new();
    let mut indices: Vec<u64> = Vec::new();
    indices.push(calculate_hash(&l_list, n as u64));
    for _ in 0..alpha {
        indices_first.push(*indices.last().unwrap() as usize);
        let mut s_ord = n.clone();
        for (i, l) in l_list.iter().enumerate() {
            let index = *indices.last().unwrap() as usize;
            let (m, p) = query_at_index(
                mtrees[i].clone(),
                mtrees[i + 1].clone(),
                points[i].clone(),
                points[i + 1].clone(),
                index,
                l_list[i],
                s_ord,
            );

            paths = [paths.clone(), m].concat();
            queried_points = [queried_points.clone(), p.clone()].concat();
            indices.push(calculate_hash(&indices, s_ord as u64));
            s_ord /= l;
        }
    }
    for _ in 0..grinding_param {
        let index = calculate_hash(&indices_first.last().unwrap(), n as u64).try_into().unwrap();
        indices_first.push(index);
        queried_points.push(points[0][index]);
        paths.push(mtrees[0].generate_proof(index).unwrap());
    }
    (paths, queried_points, indices_first)
}

pub fn fri_prove(
    f: DensePolynomial<F>, l_list: Vec<usize>, s: F, r: F, s_ord: u64, alpha: usize, grinding_param: u8,
) -> (Vec<FieldPath>, Vec<F>, Vec<Fp>, Vec<usize>) {
    let (mtrees, mroots, evals) = commit(f, l_list.clone(), s, r, s_ord.clone());

    let (paths, points, indices) = query(mtrees, evals, l_list, s_ord as usize, alpha, grinding_param);

    (paths, points, mroots, indices)
}

// Assert the equality y != \sum_i t^i*x^i*v_i
fn verify_fold_at_index(points: Vec<F>, x: F, t: F, l: usize, theta: F) -> bool {
    let y: F = points[0];
    let z_vec: Vec<F> = points[1..].to_vec();
    let mut mat: Vec<Vec<F>> = Vec::new();
    for i in 0..l {
        let mut vec_i: Vec<F> = Vec::new();
        let z: F = t.pow([i as u64]) * x;
        for j in 0..l {
            vec_i.push(z.pow([j as u64]));
        }
        mat.push(vec_i);
    }
    let g_vec = solve_linear_system(mat, z_vec);
    let mut y_supposedly: F = F::from(0);
    for (i, val) in g_vec.iter().enumerate() {
        y_supposedly += theta.pow([i as u64]) * val;
    }
    //println!("y==y_supposedly: {:?}", y==y_supposedly);
    y == y_supposedly
}

//Returns indices to query for PolyIOP together with the corresponding points
pub fn fri_verify(
    mut paths: Vec<FieldPath>, mut queried_points: Vec<F>, roots: Vec<Fp>, l_list: Vec<usize>, s: F, r: F, s_ord: u64,
    alpha: u8, grinding_param: u8,
) -> (Vec<F>, Vec<usize>) {
    let leaf_crh_params = poseidon_parameters();
    let two_to_one_params = leaf_crh_params.clone();

    let mut t_vals: Vec<F> = Vec::new();
    let mut s_vals: Vec<F> = Vec::new();
    let mut r_vals: Vec<F> = Vec::new();
    let mut s_ord_vals: Vec<u64> = Vec::new();

    s_vals.push(s);
    s_ord_vals.push(s_ord);
    r_vals.push(r);
    for (i, l) in l_list.as_slice().iter().enumerate() {
        t_vals.push(s.pow(&[s_ord / (*l as u64)]));

        r_vals.push(r_vals[i].pow(&[*l as u64]));
        s_vals.push(s_vals[i].pow(&[*l as u64]));
        s_ord_vals.push(s_ord_vals[i] / (*l as u64));
    }
    // Define all the thetas using Fiat-Shamir
    let mut theta_vec: Vec<F> = Vec::new();
    let params = poseidon_parameters();
    let mut rr: Vec<Fp> = vec![roots[0]];
    theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, rr[..1].to_vec()).unwrap(), Fp::from(0)));
    for root in roots[1..].to_vec().iter() {
        rr.push(*root);
        theta_vec.push(F::new(poseidon::CRH::<Fp>::evaluate(&params, rr.clone()).unwrap(), Fp::from(0)));
    }
    let mut indices: Vec<u64> = vec![calculate_hash(&l_list, s_ord as u64)];
    let mut i: usize = 0;
    for __ in 0..alpha {
        //points_first.push(queried_points[i]);
        for (j, l) in l_list.iter().enumerate() {
            assert!(paths[i]
                .verify(
                    &leaf_crh_params,
                    &two_to_one_params,
                    &roots[j + 1],
                    [queried_points[i].c0(), queried_points[i].c1()]
                )
                .unwrap());
            i += 1;
            for _ in 0..*l {
                assert!(paths[i]
                    .verify(
                        &leaf_crh_params,
                        &two_to_one_params,
                        &roots[j],
                        [queried_points[i].c0(), queried_points[i].c1()]
                    )
                    .unwrap());
                i += 1;
            }
        }
    }
    let mut points_first: Vec<F> = Vec::new();
    let mut indices_first: Vec<usize> = Vec::new();

    for _ in 0..alpha {
        indices_first.push(*indices.last().unwrap() as usize);
        points_first.push(queried_points[1]);
        for (i, l) in l_list.as_slice().iter().enumerate() {
            let index = *indices.last().unwrap();
            assert!(verify_fold_at_index(
                queried_points[0..*l as usize + 1].to_vec(),
                s_vals[i].pow(&[index]) * r_vals[i],
                t_vals[i],
                *l as usize,
                theta_vec[i]
            ));
            indices.push(calculate_hash(&indices, s_ord_vals[i] as u64));
            paths.drain(0..*l as usize + 1);
            queried_points.drain(0..*l as usize + 1);
        }
    }
    //queried_points = queried_points[grinding_param as usize..].to_vec();
    for i in 0..grinding_param {
        indices_first.push(calculate_hash(&indices_first.last().unwrap(), s_ord).try_into().unwrap());
        points_first.push(queried_points[i as usize]);
    }

    (points_first, indices_first)
}

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
fn calculate_hash<T: Hash>(t: &T, n: u64) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish() % n
}
