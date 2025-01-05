use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial};
use std::ops::Mul;
pub mod sidh_sike_p434;
use sidh_sike_p434::{MyFq2 as F, F as Fp};
pub mod matrix;
use matrix::*;
use merkle::{poseidon_parameters, FieldMT, FieldPath};
use std::{
    fs::File,
    io::{self, BufRead, Write},
    str::FromStr,
    time::Instant,
};

// TODO: Move to separate crate
pub mod generalized_fri;
pub mod get_roots;
pub mod isogeny_prove;
use ark_crypto_primitives::crh::poseidon;
use ark_crypto_primitives::CRHScheme;
use ark_ff::{Field, UniformRand};
use ark_std::test_rng;
use generalized_fri::*;
use isogeny_prove::{prove, verify};
pub mod merkle;
use ark_serialize::{CanonicalSerialize, Compress};

// The entry point of the application.
// Executes multiple rounds of isogeny proofs and records the results.
fn main() -> io::Result<()> {
    let mut results = Vec::new();
    // Execute 5 rounds of the isogeny proof process
    for i in 0..5 {
        if let Some(result) = round(i) {
            // Store the round index, prover time, verifier time, and proof size
            results.push((9+i, result.0, result.1, result.2));
        }
    } // Specify the path to the file where you want to save the results
    let file_path = "results_new.txt";
    write_results_to_file(&results, file_path)?;

    Ok(())
}

/// Parses the isogeny walk from a file and returns a vector of F elements.
fn parse_isogeny_walk(filename: &str) -> io::Result<Vec<F>> {
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    let mut elements = Vec::new();

    for line in reader.lines() {
        let line = line?.trim().to_string();
        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        let a: Fp;
        let b: Fp;
        if !line.contains("*a") {
            // Just a constant term
            b = Fp::from(0);
            a = Fp::from_str(&line).unwrap();
        } else if !line.contains("+") {
            // Just an 'a' term
            b = Fp::from_str(&line).unwrap();
            a = Fp::from(0);
        } else {
            // Both terms
            let (b_str, a_str) = line.split_once("*a + ").unwrap();
            a = Fp::from_str(a_str).unwrap();
            b = Fp::from_str(b_str).unwrap();
        }
        elements.push(F::new(a, b));
    }

    Ok(elements)
}
/// Writes the results of each round to a specified file.
fn write_results_to_file(results: &Vec<(usize, u64, u64, f32)>, file_path: &str) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    for result in results {
        // Write each result tuple as a comma-separated line
        writeln!(file, "{},{},{},{}", result.0, result.1, result.2, result.3)?;
    }
    Ok(())
}

/// Executes a single round of the isogeny proof process.
fn round(i: usize) -> Option<(u64, u64, f32)> {
    // We define n to be the length of the isogeny walk
    let n: usize = 2usize.pow((9 + i).try_into().unwrap());
    // Read the isogeny walk from file. Pick the first 2^{9 + i} elements from the file
    let isogeny_walk: Vec<F> = parse_isogeny_walk("isogeny_walk.txt").expect("Failed to parse file")[0..n].to_vec();
    // Define the list `l_list` based on the current round index
    // This list can contain sublists of arbitrary positive integers, as we're using a generalized
    // version of the FRI protocol which can handle arbitrary folding factors
    let l_list: Vec<usize> = vec![vec![2; 11 + i]].concat();

    // Start timing the prover
    let now = Instant::now();

    // Define the interpolation domain for the witness polynomials
    let interpolation_domain = GeneralEvaluationDomain::<F>::new(n).unwrap();
    let g = interpolation_domain.element(1);

    // Interpolate the witness polynomial as the isogeny walk over the interpolation domain
    let witness_coeffs: Vec<F> = interpolation_domain.ifft(&isogeny_walk);
    let witness: DensePolynomial<F> = DensePolynomial { coeffs: witness_coeffs };

    // Generate randomness for the blinding factor
    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);

    // Create a blinding factor polynomial
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }
        .mul(&DensePolynomial { coeffs: vec![vec![-F::from(1)], vec![F::from(0); n - 1], vec![F::from(1)]].concat() });

    // Apply the blinding factor to the witness polynomial
    let b_witness: DensePolynomial<F> = witness.clone() + blinding_factor;

    // Define the psi polynomial with evaluations (isogeny_walk[i]-isogeny_walk[i+2])^{-1} for the first n-2 elements
    // Append 2 random elements to the end of the list
    let mut psi_evals: Vec<F> =
        (0..n - 2).map(|i| (isogeny_walk[i] - isogeny_walk[i + 2]).inverse().unwrap()).collect();
    psi_evals.push(F::rand(&mut rng));
    psi_evals.push(F::rand(&mut rng));
    let psi_coeffs = interpolation_domain.ifft(&psi_evals);
    let psi: DensePolynomial<F> = DensePolynomial { coeffs: psi_coeffs };

    // Evaluate the blinded witness polynomial at specific points
    let y_start: F = b_witness.evaluate(&F::from(1));
    let y_end: F = b_witness.evaluate(&g.pow(&[n as u64 - 1]));

    // s_ord is the multiplicative order of s
    let s_ord: usize = n * 32;
    // Define the field element `r` as some element not in <s> such that r<s> is a coset
    let r: F = F::new(Fp::from(5), Fp::from(3));
    // In practice, we use <s> as the evaluation domain by transforming the polynomials P(x) -> P(r*x)
    // and then evaluate that polynomial on <s>.
    // This allows us to use FFTs for fast evaluation
    let eval_domain = GeneralEvaluationDomain::<F>::new(s_ord).unwrap();
    let s = eval_domain.element(1);

    // The below two variables are security parameteres are for the FRI protocol
    let rep_param: usize = 1;
    let grinding_param: u8 = 32;

    // Prove the validity of the isogeny walk
    let (challenge_vals, roots_fri, roots, paths_fri, points_fri, additional_paths_and_points, ws, paths_and_points) =
        prove(witness, psi, g, r, s_ord, &y_start, &y_end, l_list.clone(), rep_param, grinding_param);
    println!("Prover Time: {} s", now.elapsed().as_secs());
    let prover_time = now.elapsed().as_secs();

    let now = Instant::now();
    // Verify the proof
    let b = verify(
        challenge_vals.clone(),
        roots_fri.clone(),
        roots.clone(),
        paths_fri.clone(),
        points_fri.clone(),
        additional_paths_and_points.clone(),
        ws.clone(),
        g,
        s,
        r,
        &(n as u64),
        s_ord,
        &y_start,
        &y_end,
        l_list,
        rep_param,
        grinding_param,
        paths_and_points.clone(),
    );
    println!("Verifier Time: {} ms", now.elapsed().as_millis());
    let verifier_time = now.elapsed().as_millis() as u64;
    // Compute the total size of all the proof elements
    if b {
        let size1 = challenge_vals.serialized_size(Compress::Yes);
        let size2 = roots.serialized_size(Compress::Yes);
        let size3 = points_fri.serialized_size(Compress::Yes);
        let size4 = roots_fri.serialized_size(Compress::Yes);
        let size5 = paths_fri.serialized_size(Compress::Yes);
        let size6 = additional_paths_and_points[0].paths_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[0].points_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[0].paths_minus.serialized_size(Compress::Yes)
            + additional_paths_and_points[0].points_minus.serialized_size(Compress::Yes);
        let size7 = additional_paths_and_points[1].paths_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[1].points_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[1].paths_minus.serialized_size(Compress::Yes)
            + additional_paths_and_points[1].points_minus.serialized_size(Compress::Yes);
        let size8 = additional_paths_and_points[2].paths_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[2].points_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[2].paths_minus.serialized_size(Compress::Yes)
            + additional_paths_and_points[2].points_minus.serialized_size(Compress::Yes);
        let size9 = paths_and_points[0].root_f.serialized_size(Compress::Yes)
            + paths_and_points[0].root_g.serialized_size(Compress::Yes)
            + paths_and_points[0].root_u.serialized_size(Compress::Yes)
            + paths_and_points[0].merkle_paths_f.serialized_size(Compress::Yes)
            + paths_and_points[0].merkle_paths_g.serialized_size(Compress::Yes)
            + paths_and_points[0].merkle_paths_u.serialized_size(Compress::Yes)
            + paths_and_points[0].queried_points_f.serialized_size(Compress::Yes)
            + paths_and_points[0].queried_points_u.serialized_size(Compress::Yes)
            + paths_and_points[0].queried_points_g.serialized_size(Compress::Yes);
        let size10 = paths_and_points[1].root_f.serialized_size(Compress::Yes)
            + paths_and_points[1].root_g.serialized_size(Compress::Yes)
            + paths_and_points[1].root_u.serialized_size(Compress::Yes)
            + paths_and_points[1].merkle_paths_f.serialized_size(Compress::Yes)
            + paths_and_points[1].merkle_paths_g.serialized_size(Compress::Yes)
            + paths_and_points[1].merkle_paths_u.serialized_size(Compress::Yes)
            + paths_and_points[1].queried_points_f.serialized_size(Compress::Yes)
            + paths_and_points[1].queried_points_u.serialized_size(Compress::Yes)
            + paths_and_points[1].queried_points_g.serialized_size(Compress::Yes);
        let size11 = paths_and_points[2].root_f.serialized_size(Compress::Yes)
            + paths_and_points[2].root_g.serialized_size(Compress::Yes)
            + paths_and_points[2].root_u.serialized_size(Compress::Yes)
            + paths_and_points[2].merkle_paths_f.serialized_size(Compress::Yes)
            + paths_and_points[2].merkle_paths_g.serialized_size(Compress::Yes)
            + paths_and_points[2].merkle_paths_u.serialized_size(Compress::Yes)
            + paths_and_points[2].queried_points_f.serialized_size(Compress::Yes)
            + paths_and_points[2].queried_points_u.serialized_size(Compress::Yes)
            + paths_and_points[2].queried_points_g.serialized_size(Compress::Yes);

        println!(
            "Proof Size: {} kB",
            ((size1 + size2 + size3 + size4 + size5 + size6 + size7 + size8 + size9 + size10 + size11) as f32)
                / 1000f32
        );
        let proof_size = ((size1 + size2 + size3 + size4 + size5 + size6 + size7 + size8 + size9 + size10 + size11)
            as f32)
            / 1000f32;
        println!("Verification successful");
        Some((prover_time, verifier_time, proof_size))
    } else {
        println!("Verification failed");
        None
    }
}
