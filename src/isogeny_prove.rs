use std::ops::{Div, Sub, Mul};

use super::*;
use ark_ff::UniformRand;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_std::{rand::Rng, test_rng};
use merkle::{poseidon_parameters, FieldMT};

// Commit to a randomized vector based on polynomial evaluations in order to preserve the 
// SHVZK property of the FRI protocol
fn commit_with_rand(f: DensePolynomial<F>, eval_domain: Vec<F>) -> (Vec<F>, FieldMT, Vec<F>) {
    let s_ord: u64 = eval_domain.len() as u64;
    let params = poseidon_parameters();
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();

    // Generate a random vector of F elements
    let mut rng = test_rng();
    let random_vec: Vec<F> = (0..s_ord).into_iter().map(|_| F::rand(&mut rng)).collect();
    let pure_evals: Vec<F> = eval_domain.iter().map(|x| f.evaluate(x)).collect();
    // Generate evaluation vector which adds f(eval_vec[i]) and random_vec[i]
    let randomized_evals: Vec<F> = pure_evals.iter().zip(random_vec.iter()).map(|(x, r)| x + r).collect();

    let evals_slice: Vec<Vec<Fp>> = randomized_evals.iter().map(|x| vec![x.c0(), x.c1()]).collect();
    let random_vec_slice: Vec<Vec<Fp>> = random_vec.iter().map(|x| vec![x.c0(), x.c1()]).collect();
    
    let randomized_evals = vec![
        evals_slice,
        random_vec_slice,
    ]
    .concat();

    (pure_evals, FieldMT::new(&leaf_crh_params, &two_to_one_params, randomized_evals).unwrap(), random_vec)
}
// Function to compute pseudorandom list of n indices in span [0, 2^k] given some input c
fn compute_pseudorandom_indices(c: F, n: usize, k: u32) -> Vec<usize> {
    let mut indices = Vec::with_capacity(n);
    let _current_hash = c;

    for _ in 0..n {
        // todo: Make pseudorandom
        // Generate random index in [0, 2^k)
        let mut rng = test_rng();
        let index: usize = rng.gen_range(0..2u64.pow(k)).try_into().unwrap();
        indices.push(index);
    }

    indices
}

// Executes steps 1.-6. of the BlindedEval subprotocol of FRI
fn blinded_eval_first_six_steps(
    evals_vec: Vec<F>, merkle_tree: FieldMT, random_vec: Vec<F>, eval_domain: Vec<F>, deg: usize, z: F, beta: usize,
) -> (F, PathsAndPoints, DensePolynomial<F>, FieldMT, Vec<F>) {
    // First define random polynomial g with degree deg
    let mut rng = test_rng();
    let g_coeffs: Vec<F> = (0..deg + 1).into_iter().map(|_| F::rand(&mut rng)).collect();
    let g: DensePolynomial<F> = DensePolynomial { coeffs: g_coeffs };
    let params = poseidon_parameters();
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();

    let s_ord = eval_domain.len();
    let k: u32 = (((s_ord as f32).log2()).ceil()) as u32;

    // Commit to g with rand
    let (g_evals, g_mtree, random_vec_g) = commit_with_rand(g.clone(), eval_domain.clone());
    let g_eval = g.evaluate(&z);

    // Generate random challenge from merkle root
    let c: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![g_mtree.root()]).unwrap(), Fp::from(0));
    // Define H as g_evals + c*evals_vec
    let H: Vec<F> = g_evals.iter().zip(evals_vec.iter()).map(|(g_eval, eval)| *g_eval + c * *eval).collect();

    // Define T as random_vec_g + c*random_vec
    let T: Vec<F> = random_vec_g.iter().zip(random_vec.iter()).map(|(rand_g, rand)| *rand_g + c * *rand).collect();

    // Define U as H+T concatenated with T
    let mut U: Vec<F> = H.iter().zip(T.iter()).map(|(h, t)| *h + *t).collect();
    U.extend(T.iter().cloned());
    let u_slice: Vec<Vec<Fp>> = U.iter().map(|x| vec![x.c0(), x.c1()]).collect();
    // Commit to U
    let u_mtree = FieldMT::new(&leaf_crh_params, &two_to_one_params, u_slice).unwrap();

    // Compute list of beta indices based on c
    let beta_indices = compute_pseudorandom_indices(c, beta, k);

    // Define vecs of merkle paths for merkle_tree, merkle_tree_u, and merkle_tree_g
    let mut merkle_paths: Vec<FieldPath> = vec![];
    let mut merkle_paths_u: Vec<FieldPath> = vec![];
    let mut merkle_paths_g: Vec<FieldPath> = vec![];

    // Define vecs of queried points
    let mut queried_points: Vec<F> = vec![];
    let mut queried_points_u: Vec<F> = vec![];
    let mut queried_points_g: Vec<F> = vec![];
    // For each index, query the merkle trees merkle_tree, merkle_tree_u and merkle_tree_g
    for index in beta_indices.iter() {
        let u_eval = U[*index];
        let g_eval = g_evals[*index] + random_vec_g[*index];
        let f_eval = evals_vec[*index] + random_vec[*index];
        queried_points.push(f_eval);
        queried_points_u.push(u_eval);
        queried_points_g.push(g_eval);

        merkle_paths.push(merkle_tree.generate_proof(*index).unwrap());
        merkle_paths_u.push(u_mtree.generate_proof(*index).unwrap());
        merkle_paths_g.push(g_mtree.generate_proof(*index).unwrap());
    }

    // Put all relevant proof elements into a PathsAndPoints data structure
    let paths_and_points = PathsAndPoints {
        root_f: merkle_tree.root(),
        root_g: g_mtree.root(),
        root_u: u_mtree.root(),

        merkle_paths_f: merkle_paths,
        merkle_paths_u: merkle_paths_u,
        merkle_paths_g: merkle_paths_g,

        queried_points_f: queried_points,
        queried_points_u: queried_points_u,
        queried_points_g: queried_points_g,
    };
    (g_eval, paths_and_points, g, u_mtree, U)
}

// Verifies steps 1.-6. of the BlindedEval subprotocol of FRI
fn verify_blinded_eval_first_six_steps(_g_eval: F, paths_and_points: &PathsAndPoints, _z: F, beta: usize) -> bool {
    let mut result: bool = true;

    // Define the parameters for hashing
    let params = poseidon_parameters();
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();

    let root_f = paths_and_points.root_f;
    let root_g = paths_and_points.root_g;
    let c: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![root_g]).unwrap(), Fp::from(0));

    let root_u = paths_and_points.root_u;
    for i in 0..beta {
        let _path_f = &paths_and_points.merkle_paths_f[i];
        let point_f = &paths_and_points.queried_points_f[i];

        let _path_g = &paths_and_points.merkle_paths_g[i];
        let point_g = &paths_and_points.queried_points_g[i];

        let _path_u = &paths_and_points.merkle_paths_u[i];
        let point_u = &paths_and_points.queried_points_u[i];

        // Verify that all the Merkle paths are valid and the relation U == G + c * F holds at the queried points
        result = result
            && paths_and_points.merkle_paths_f[i]
                .verify(&leaf_crh_params, &two_to_one_params, &root_f, [point_f.c0(), point_f.c1()])
                .unwrap()
            && paths_and_points.merkle_paths_g[i]
                .verify(&leaf_crh_params, &two_to_one_params, &root_g, [point_g.c0(), point_g.c1()])
                .unwrap()
            && paths_and_points.merkle_paths_u[i]
                .verify(&leaf_crh_params, &two_to_one_params, &root_u, [point_u.c0(), point_u.c1()])
                .unwrap()
            && (*point_u == (point_g + &(c * point_f)));
    }
    result
}

// Define struct for three vecs of merkle paths and three vecs of queried points
#[derive(Clone)]
pub struct PathsAndPoints {
    pub root_f: Fp,
    pub root_g: Fp,
    pub root_u: Fp,

    pub merkle_paths_f: Vec<FieldPath>,
    pub merkle_paths_u: Vec<FieldPath>,
    pub merkle_paths_g: Vec<FieldPath>,

    pub queried_points_f: Vec<F>,
    pub queried_points_u: Vec<F>,
    pub queried_points_g: Vec<F>,
}
// Define struct for paths and points used as part of FRI
#[derive(Clone)]
pub struct PathsPointsPlusMinus {
    pub root: Fp,

    pub paths_plus: Vec<FieldPath>,
    pub points_plus: Vec<F>,
    pub paths_minus: Vec<FieldPath>,
    pub points_minus: Vec<F>,
}
// This function finds the points and corresponding Merkle paths of the indices
fn query_relevant_indices(
    points_vec: Vec<F>, merkle_tree: FieldMT, indices: Vec<usize>, n: usize,
) -> PathsPointsPlusMinus {
    let mut queried_points_plus: Vec<F> = vec![];
    let mut queried_paths_plus: Vec<FieldPath> = vec![];
    let mut queried_points_minus: Vec<F> = vec![];
    let mut queried_paths_minus: Vec<FieldPath> = vec![];
    for i in indices {
        queried_points_minus.push(points_vec[i]);
        queried_paths_minus.push(merkle_tree.generate_proof(i).unwrap());

        let i_plus = i + n;
        queried_points_plus.push(points_vec[i_plus]);
        queried_paths_plus.push(merkle_tree.generate_proof(i_plus).unwrap());
    }
    PathsPointsPlusMinus {
        root: merkle_tree.root(),
        paths_plus: queried_paths_plus,
        points_plus: queried_points_plus,
        paths_minus: queried_paths_minus,
        points_minus: queried_points_minus,
    }
}

// Verifies the Merkle paths, outputs the value m - m_{+} as per step 7. of FRI.BatchedBlindedEval
fn verify_paths_and_points(paths_and_points: Vec<PathsPointsPlusMinus>) -> Option<Vec<Vec<F>>> {
    let params = poseidon_parameters();
    let leaf_crh_params = params.clone();
    let two_to_one_params = params.clone();

    let mut result: Vec<Vec<F>> = vec![];
    for i in 0..paths_and_points.len() {
        let current_pp = &paths_and_points[i];
        let root = current_pp.root;
        let points_plus = &current_pp.points_plus;
        let points_minus = &current_pp.points_minus;
        let paths_plus = &current_pp.paths_plus;
        let paths_minus = &current_pp.paths_minus;
        let mut temp_result: Vec<F> = vec![];
        for j in 0..points_plus.len() {
            let point_plus: F = points_plus[j];
            let point_minus: F = points_minus[j];
            let b_plus = paths_plus[j]
                .verify(&leaf_crh_params, &two_to_one_params, &root, [point_plus.c0(), point_plus.c1()])
                .unwrap();
            let b_minus = paths_minus[j]
                .verify(&leaf_crh_params, &two_to_one_params, &root, [point_minus.c0(), point_minus.c1()])
                .unwrap();
            // If either of the paths is invalid, exit function and return None
            if !b_plus || !b_minus {
                return None;
            }
            // If the Merkle paths are both correct, push the difference between the points to the
            // output vector
            temp_result.push(point_minus - point_plus);
        }
        result.push(temp_result);
    }
    Some(result)
}
// Witness is the witness polynomial, psi the inverse of w(x)-w(g^{2}*x), g the generator of the interpolation domain,
//the evaluation domain is r<s>. Finally, s_ord is the size of E. 
pub fn prove(
    witness: DensePolynomial<F>, psi: DensePolynomial<F>, g: F, s: F, r: F, s_ord: u64, y_start: &F, y_end: &F,
    l_list: Vec<usize>, rep_param: usize, grinding_param: u8,
) -> (Vec<F>, Vec<Fp>, Vec<Fp>, Vec<FieldPath>, Vec<F>, Vec<PathsPointsPlusMinus>, Vec<F>, Vec<PathsAndPoints>) {
    let n: usize = witness.coeffs.len();

    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    // Blind the witness
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }.mul(&DensePolynomial {
        coeffs: vec![vec![-F::from(1)], vec![F::from(0); n - 1], vec![F::from(1)]].concat(),
    });

    let b_witness: DensePolynomial<F> = witness.clone() + blinding_factor.clone();
    let b_witness_plus: DensePolynomial<F> = DensePolynomial {
        coeffs: b_witness.coeffs.iter().enumerate().map(|(i, coeff)| coeff * &g.pow(&[i as u64])).collect(),
    };
    let b_witness_plus_plus: DensePolynomial<F> = DensePolynomial {
        coeffs: b_witness_plus.coeffs.iter().enumerate().map(|(i, coeff)| coeff * &g.pow(&[i as u64])).collect(),
    };
    // Commit to the blinded witness and psi
    let params = poseidon_parameters();
    let _leaf_crh_params = params.clone();
    let _two_to_one_params = params.clone();

    // Define the evaluation domain 
    let eval_domain: Vec<F> = (0..s_ord).into_iter().map(|i| r * s.pow([i])).collect();

    // Commit with randomness to the blinded witness and psi
    let (witness_evals, witness_mtree, random_vec_witness) = commit_with_rand(b_witness.clone(), eval_domain.clone());
    let (psi_evals, psi_mtree, random_vec_psi) = commit_with_rand(psi.clone(), eval_domain.clone());

    // Store the Merkle roots in a vector
    let mut roots: Vec<Fp> = vec![witness_mtree.root(), psi_mtree.root()];

    // Use the roots as a pseudorandom generator for the \alpha_{i} used in the isogeny walk protocol
    let alpha_1: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots[..2].to_vec().clone()).unwrap(), Fp::from(0));
    let alpha_2: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_1.c0()]).unwrap(), Fp::from(0));
    let alpha_3: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_2.c0()]).unwrap(), Fp::from(0));
    let alpha_4: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_3.c0()]).unwrap(), Fp::from(0));

    // Compute C_{i}(x)
    let c1: DensePolynomial<F> = initial_poly(&y_start, b_witness.clone());
    let c2: DensePolynomial<F> = mod_poly_poly(&b_witness, &b_witness_plus, n, g);
    let c3: DensePolynomial<F> = psi_poly(&b_witness, &b_witness_plus_plus, &psi, n, g);
    let c4: DensePolynomial<F> = final_poly(&y_end, b_witness.clone(), g, n as u64);

    // Compute C(x)
    let c: DensePolynomial<F> = compute_c(c1, c2, c3, c4, &vec![alpha_1, alpha_2, alpha_3, alpha_4], &n);
    
    // Commit to C(x)
    let (c_evals, c_mtree, random_vec_c) = commit_with_rand(c.clone(), eval_domain.clone());

    roots.push(c_mtree.root());

    // Compute the pseudorandom evaluation point z
    let z: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0));

    let gz: F = g * z;
    let ggz: F = g * gz;

    // Evaluate the polynomials on z 
    let witness_y: F = b_witness.evaluate(&z);
    let witness_y_plus: F = b_witness.evaluate(&gz);
    let witness_y_plus_plus: F = b_witness.evaluate(&ggz);
    let psi_y: F = psi.evaluate(&z);
    let c_y: F = c.evaluate(&z);

    // Define E as in the batched general FRI protocol
    let E: usize = 4 * n;

    // Execute the first 6 steps of the generalized FRI BlindedEval protocol on \phi(x), \psi(x), C(x)
    let (witness_g_eval, witness_paths_and_points, witness_g, witness_u_mtree, witness_u) =
        blinded_eval_first_six_steps(
            witness_evals,
            witness_mtree,
            random_vec_witness.clone(),
            eval_domain.clone(),
            n,
            z,
            rep_param,
        );
    let (psi_g_eval, psi_paths_and_points, psi_g, psi_u_mtree, psi_u) = blinded_eval_first_six_steps(
        psi_evals,
        psi_mtree,
        random_vec_psi.clone(),
        eval_domain.clone(),
        n,
        z,
        rep_param,
    );
    let (c_g_eval, c_paths_and_points, c_g, c_u_mtree, c_u) =
        blinded_eval_first_six_steps(c_evals, c_mtree, random_vec_c.clone(), eval_domain.clone(), n, z, rep_param);
    let challenge_vals: Vec<F> = vec![witness_y, witness_y_plus, witness_y_plus_plus, psi_y, c_y];

    // Save the evals vectors to a vector
    let ws = vec![witness_g_eval, psi_g_eval, c_g_eval];

    // Finally, define and commit to the composition polynomial P(x)
    let mut zeta_vec: Vec<Fp> = vec![poseidon::CRH::<Fp>::evaluate(&params, vec![z.c0()]).unwrap()];
    for _ in 0..4 {
        zeta_vec.push(poseidon::CRH::<Fp>::evaluate(&params, zeta_vec.clone()).unwrap());
    }
    let const_witness =
        F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![witness_paths_and_points.root_g]).unwrap(), Fp::from(0));
    let const_psi =
        F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![psi_paths_and_points.root_g]).unwrap(), Fp::from(0));
    let const_c = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![c_paths_and_points.root_g]).unwrap(), Fp::from(0));

    let Q_0 =
        DensePolynomial { coeffs: vec![vec![F::from(0); E - n], vec![F::new(zeta_vec[0], Fp::from(0))]].concat() }
            .mul(
                &(witness_g.clone()
                    + DensePolynomial { coeffs: vec![const_witness] }.mul(&b_witness)
                    + DensePolynomial { coeffs: vec![-ws[0] - const_witness * challenge_vals[0]] }),
            )
            .div(&DensePolynomial { coeffs: vec![-z, F::from(1)] });
    let _Q_1 =
    DensePolynomial { coeffs: vec![vec![F::from(0); E - n], vec![F::new(zeta_vec[1], Fp::from(0))]].concat() }
        .mul(
            &(witness_g.clone()
                + DensePolynomial { coeffs: vec![const_witness] }.mul(&b_witness)
                + DensePolynomial { coeffs: vec![-ws[0] - const_witness * challenge_vals[1]] }),
        )
        .div(&DensePolynomial { coeffs: vec![-gz, F::from(1)] });
        let _Q_2 =
        DensePolynomial { coeffs: vec![vec![F::from(0); E - n], vec![F::new(zeta_vec[2], Fp::from(0))]].concat() }
            .mul(
                &(witness_g.clone()
                    + DensePolynomial { coeffs: vec![const_witness] }.mul(&b_witness)
                    + DensePolynomial { coeffs: vec![-ws[0] - const_witness * challenge_vals[2]] }),
            )
            .div(&DensePolynomial { coeffs: vec![-ggz, F::from(1)] });
    let Q_3 =
        DensePolynomial { coeffs: vec![vec![F::from(0); E - n - 1], vec![F::new(zeta_vec[3], Fp::from(0))]].concat() }
            .mul(
                &(psi_g.clone()
                    + DensePolynomial { coeffs: vec![const_psi] }.mul(&psi)
                    + DensePolynomial { coeffs: vec![-ws[1] - const_psi * challenge_vals[3]] }),
            )
            .div(&DensePolynomial { coeffs: vec![-z, F::from(1)] });
    let Q_4 = (c_g.clone()
        + DensePolynomial { coeffs: vec![const_c] }.naive_mul(&c)
        + DensePolynomial { coeffs: vec![-ws[2] - const_c * challenge_vals[4]] })
    .div(&DensePolynomial { coeffs: vec![-z, F::from(1)] });
    
    let composition_poly = 
        Q_0 
        //+ Q_1
        //+ Q_2
        + Q_3
        + Q_4;

    // Perform the FRI protocol on the composition polynomial P(x)
    let (paths_fri, points_fri, roots_fri, indices) =
        fri_prove(composition_poly, l_list, s, r, s_ord, rep_param, grinding_param);
    let two_k = random_vec_witness.clone().len();

    let witness_paths_points_plus_minus =
        query_relevant_indices(witness_u.clone(), witness_u_mtree.clone(), indices.clone(), two_k);
    let psi_paths_points_plus_minus = query_relevant_indices(psi_u, psi_u_mtree, indices.clone(), two_k);
    let c_paths_points_plus_minus = query_relevant_indices(c_u, c_u_mtree, indices.clone(), two_k);
    let paths_and_points = vec![witness_paths_and_points, psi_paths_and_points, c_paths_and_points];
    let additional_paths_and_points: Vec<PathsPointsPlusMinus> =
        vec![witness_paths_points_plus_minus, psi_paths_points_plus_minus, c_paths_points_plus_minus];

    (challenge_vals, roots_fri, roots, paths_fri, points_fri, additional_paths_and_points, ws, paths_and_points)
}

// Verifies a proof for the l-isogeny walk protocol
pub fn verify(
    challenges: Vec<F>, roots_fri: Vec<Fp>, roots: Vec<Fp>, paths_fri: Vec<FieldPath>, points_fri: Vec<F>,
    additional_paths_and_points: Vec<PathsPointsPlusMinus>, ws: Vec<F>, g: F, s: F, r: F, n: &u64, s_ord: u64,
    y_start: &F, y_end: &F, l_list: Vec<usize>, rep_param: usize, grinding_param: u8,
    paths_and_points: Vec<PathsAndPoints>,
) -> bool {
    // Compute z, alphas and zetas
    let params = poseidon_parameters();
    let _leaf_crh_params = params.clone();
    let _two_to_one_params = params.clone();
    let z: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots.clone()).unwrap(), Fp::from(0));
    let gz: F = g * z;
    let _ggz: F = g * gz;
    let alpha_1: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, roots[..2].to_vec().clone()).unwrap(), Fp::from(0));
    let alpha_2: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_1.c0()]).unwrap(), Fp::from(0));
    let alpha_3: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_2.c0()]).unwrap(), Fp::from(0));
    let alpha_4: F = F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![alpha_3.c0()]).unwrap(), Fp::from(0));
    let mut zeta_vec: Vec<Fp> = vec![poseidon::CRH::<Fp>::evaluate(&params, vec![z.c0()]).unwrap()];

    for _ in 0..4 {
        zeta_vec.push(poseidon::CRH::<Fp>::evaluate(&params, zeta_vec.clone()).unwrap());
    }
    // Check that the FRI queries are correct
    let (points_first, indices_first) = fri_verify(
        paths_fri,
        points_fri.clone(),
        roots_fri,
        l_list.clone(),
        s.clone(),
        r.clone(),
        s_ord.clone(),
        rep_param as u8,
        grinding_param,
    );
    // Check that the challenges were computed correctly
    let E: u64 = 4 * n;
    let c1: F = initial_challenge(y_start, &challenges[0], &z);
    let c2: F = mod_challenge(&challenges[0], &challenges[1], &z, &g, &n);
    let c3: F = psi_challenge(&challenges[0], &challenges[2], &challenges[3], &z, n, &g);
    let c4: F = final_challenge(y_end, &challenges[0], &z, n, &g);

    let asserted_c: F = alpha_1 * z.pow(&[E - n - 2]) * c1
        + alpha_2 * z.pow(&[E - 3 * n - 13]) * c2
        + alpha_3 * z.pow(&[E - n - 5]) * c3
        + alpha_4 * z.pow(&[E - n - 2]) * c4;
    assert_eq!(asserted_c, challenges[4]);

    // Check consistency between P(x) in FRI and the committed-to polynomials
    let const_witness =
        F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![paths_and_points.clone()[0].root_g]).unwrap(), Fp::from(0));
    let const_psi =
        F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![paths_and_points.clone()[1].root_g]).unwrap(), Fp::from(0));
    let const_c =
        F::new(poseidon::CRH::<Fp>::evaluate(&params, vec![paths_and_points.clone()[2].root_g]).unwrap(), Fp::from(0));

    let b_0 = verify_blinded_eval_first_six_steps(ws[0], &paths_and_points[0], z, rep_param);
    let b_1 = verify_blinded_eval_first_six_steps(ws[1], &paths_and_points[1], z, rep_param);
    let b_2 = verify_blinded_eval_first_six_steps(ws[2], &paths_and_points[2], z, rep_param);
    if !(b_0 && b_1 && b_2) {
        return false;
    }

    // If the Merkle paths correctly verify, check the polynomial relation
    if let Some(points_for_consistency) = verify_paths_and_points(additional_paths_and_points) {
        let numerator_0 = ws[0] + const_witness * challenges[0];
        //let numerator_1 = ws[0] + const_witness * challenges[1];
        //let numerator_2 = ws[0] + const_witness * challenges[2];
        let numerator_3 = ws[1] + const_psi * challenges[3];
        let numerator_4 = ws[2] + const_c * challenges[4];

        for (i, index) in indices_first.iter().enumerate() {
            let x_0: F = r * s.pow(&[*index as u64]);
            let witness_val_1: F = (points_for_consistency[0][i] - numerator_0) / (x_0 - z);
            //let witness_val_2: F = (points_for_consistency[0][i] - numerator_1) / (x_0 - gz);
            //let witness_val_3: F = (points_for_consistency[0][i] - numerator_2) / (x_0 - ggz);
            let psi_val: F = (points_for_consistency[1][i] - numerator_3) / (x_0 - z);
            let c_val: F = (points_for_consistency[2][i] - numerator_4) / (x_0 - z);
            let asserted_p: F = F::new(zeta_vec[0], Fp::from(0)) * x_0.pow(&[E - n]) * witness_val_1
                //+ F::new(zeta_vec[1], Fp::from(0)) * x_0.pow(&[E - n - 1]) * witness_val_2
                //+ F::new(zeta_vec[2], Fp::from(0)) * x_0.pow(&[E - n - 1]) * witness_val_3
                + F::new(zeta_vec[3], Fp::from(0)) * x_0.pow(&[E - n-1]) * psi_val
                + c_val;
            if asserted_p != points_first[i] {
                return false;
            }
        }
        return true;
    }
    false
}

// The challenge polynomial based on the second modular polynomial
fn mod_challenge(x: &F, y: &F, z: &F, g: &F, T: &u64) -> F {
    let eval: F = x * x * x + y * y * y - x * x * y * y + F::from(1488u128) * (x * x * y + y * y * x)
        - F::from(162000u128) * (x * x + y * y)
        + F::from(40773375u128) * x * y
        + F::from(8748000000u128) * (x + y)
        - F::from(157464000000000u128);

    (*&z - &g.pow(&[*T - 1])) * eval / (*&z.pow(&[*T]) - &F::from(1))
}

// Returns (p(x)-y_0)/(x - 1)
fn initial_poly(y_0: &F, p: DensePolynomial<F>) -> DensePolynomial<F> {
    (p + DensePolynomial { coeffs: vec![-*y_0] }).div(&DensePolynomial { coeffs: vec![F::from(-1), F::from(1)] })
}
//Returns (p(x)-y_end)/(x - g^(T-1))
fn final_poly(y_end: &F, p: DensePolynomial<F>, g: F, T: u64) -> DensePolynomial<F> {
    (p + DensePolynomial { coeffs: vec![-*y_end] }).div(&DensePolynomial { coeffs: vec![-g.pow(&[T - 1]), F::from(1)] })
}

// Returns the composite polynomial
fn compute_c(
    c1: DensePolynomial<F>, c2: DensePolynomial<F>, c3: DensePolynomial<F>, c4: DensePolynomial<F>, alphas: &Vec<F>,
    T: &usize,
) -> DensePolynomial<F> {
    let E: usize = 4 * T;
    let deg_1: usize = E - T - 2;
    let deg_2: usize = E - T - 5;
    let deg_3: usize = E - 3 * T - 13;

    c1.mul(&DensePolynomial { coeffs: vec![vec![F::from(0); deg_1], vec![alphas[0]]].concat() })
        + c2.mul(&DensePolynomial { coeffs: vec![vec![F::from(0); deg_3], vec![alphas[1]]].concat() })
        + c3.mul(&DensePolynomial { coeffs: vec![vec![F::from(0); deg_2], vec![alphas[2]]].concat() })
        + c4.mul(&DensePolynomial { coeffs: vec![vec![F::from(0); deg_1], vec![alphas[3]]].concat() })
}
// Returns (x-g^(T-1))*Phi_2(p(x), q(x))/(x^n - 1)
fn mod_poly_poly(p: &DensePolynomial<F>, q: &DensePolynomial<F>, T: usize, g: F) -> DensePolynomial<F> {
    let p_squared: DensePolynomial<F> = p.mul(p);
    let q_squared: DensePolynomial<F> = q.mul(q);
    let p_cubed: DensePolynomial<F> = p_squared.mul(p);
    let q_cubed: DensePolynomial<F> = q_squared.mul(q);
    let p_squared_q: DensePolynomial<F> = p_squared.mul(q);
    let q_squared_p: DensePolynomial<F> = q_squared.mul(p);
    let pq: DensePolynomial<F> = p.mul(q);

    let temp: DensePolynomial<F> = p_cubed
        + q_cubed.sub(&p_squared_q.mul(q))
        + DensePolynomial { coeffs: vec![F::from(1488u128)] }.mul(&(p_squared_q + q_squared_p))
        + DensePolynomial { coeffs: vec![-F::from(162000u128)] }.mul(&(p_squared + q_squared))
        + DensePolynomial { coeffs: vec![F::from(40773375u128)] }.mul(&pq)
        + DensePolynomial { coeffs: vec![F::from(8748000000u128)] }.mul(&(p + q))
        + DensePolynomial { coeffs: vec![-F::from(157464000000000u128)] };

    temp.mul(&DensePolynomial { coeffs: vec![-g.pow(&[T as u64 - 1]), F::from(1)] })
        .div(&DensePolynomial { coeffs: [vec![-F::from(1)], vec![F::from(0); T - 1], vec![F::from(1)]].concat() })
}

fn initial_challenge(y_0: &F, eval: &F, x_0: &F) -> F {
    (eval - y_0) / (x_0 - &F::from(1))
}
fn final_challenge(y_end: &F, eval: &F, x_0: &F, n: &u64, g: &F) -> F {
    (eval - y_end) / (x_0 - &g.pow(&[*n - 1]))
}
//Returns ((x-g^(T-2))*(x-g^(T-1))*(p(x)-q(x))*psi(x)-1)/(x^T-1)
fn psi_poly(
    p: &DensePolynomial<F>, q: &DensePolynomial<F>, psi: &DensePolynomial<F>, T: usize, g: F,
) -> DensePolynomial<F> {
    let diff: DensePolynomial<F> = p.sub(q);
    let g_pow: F = g.pow(&[T as u64 - 2]);
    let x_1_poly: DensePolynomial<F> = DensePolynomial { coeffs: vec![-g_pow, F::from(1)] }
        .mul(&DensePolynomial { coeffs: vec![-g_pow * g, F::from(1)] });

    x_1_poly
        .mul(&(diff.mul(psi) + DensePolynomial { coeffs: vec![F::from(-1)] }))
        .div(&DensePolynomial { coeffs: [vec![-F::from(1)], vec![F::from(0); T - 1], vec![F::from(1)]].concat() })
}

fn psi_challenge(y_witness: &F, y_witness_plusplus: &F, y_psi: &F, x_0: &F, n: &u64, g: &F) -> F {
    let g_pow: F = g.pow(&[*n as u64 - 2]);
    let g_prefactor: F = (*&x_0 - &g_pow) * (*&x_0 - &(g_pow * g));

    g_prefactor * ((y_witness - y_witness_plusplus) * y_psi - &F::from(1)) / (*&x_0.pow(&[*n]) - &F::from(1))
}
