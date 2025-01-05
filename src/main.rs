use ark_ff::MontFp;
use ark_poly::univariate::DensePolynomial;
//pub mod field;
pub mod sidh_sike_p434;
use sidh_sike_p434::{Fq2 as F, F as Fp};
pub mod matrix;
use matrix::*;
use merkle::{poseidon_parameters, FieldMT, FieldPath};
use std::{
    fs::{self, File},
    io::{self, Write},
    path::Path,
    str::FromStr,
    time::Instant,
};

// TODO: Move to separate crate
pub mod generalized_fri;
pub mod get_roots;
pub mod isogeny_prove;
use ark_crypto_primitives::crh::poseidon;
use ark_crypto_primitives::CRHScheme;
use ark_ff::Field;
use ark_ff::UniformRand;
use ark_poly::Polynomial;
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
            results.push((i, result.0, result.1, result.2));
        }
    } // Specify the path to the file where you want to save the results
    let file_path = "results_new.txt";
    write_results_to_file(&results, file_path)?;

    Ok(())
}

/// Reads lines from a file and parses them into a vector of field elements.
///
/// # Arguments
///
/// * `filename` - The path to the file to read.
///
/// # Returns
///
/// A `Result` containing a vector of `F` elements or an `io::Error`.
fn lines_from_file(filename: impl AsRef<Path>) -> io::Result<Vec<F>> {
    let file_content = fs::read_to_string(filename)?;
    let lines = file_content.split(',').collect::<Vec<_>>();
    let mut results = Vec::new();
    for line in lines {
        let a: Fp;
        let b: Fp;
        // If the line does not contain "*a", parse it as a single Fp element
        if !line.contains("*a") {
            println!("line: {:?}", line);
            a = Fp::from_str(line).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp4"))?;
            b = Fp::from(0);
        } else if !line.contains("+") {
            // If the line contains "*a" but not "+", parse accordingly
            let mut parts = line.trim().split("*a");
            b = Fp::from_str(parts.next().unwrap().trim())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp2"))?;
            a = Fp::from(0);
        } else {
            // If the line contains both "*a" and "+", parse both parts
            let mut parts = line.trim().split("*a +");
            b = Fp::from_str(parts.next().unwrap())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp1"))?;
            a = Fp::from_str(parts.next().unwrap().trim())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp5"))?;
        }
        // Combine the parsed parts into a quadratic extension field element
        results.push(F::new(a, b));
    }
    Ok(results)
}

/// Writes the results of each round to a specified file.
///
/// # Arguments
///
/// * `results` - Reference to a vector of tuples containing round results.
/// * `file_path` - The path to the file where results will be written.
///
/// # Returns
///
/// A `Result` indicating success or an `io::Error`.
fn write_results_to_file(results: &Vec<(usize, u64, u64, f32)>, file_path: &str) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    for result in results {
        // Write each result tuple as a comma-separated line
        writeln!(file, "{},{},{},{}", result.0, result.1, result.2, result.3)?;
    }
    Ok(())
}

/// Executes a single round of the isogeny proof process.
///
/// # Arguments
///
/// * `i` - The index of the current round.
///
/// # Returns
///
/// An `Option` containing a tuple of prover time, verifier time, and proof size if verification succeeds.
fn round(i: usize) -> Option<(u64, u64, f32)> {
    // Start timing the prover
    let now = Instant::now();

    // Define the list `l_list` based on the current round index
    // This list can contain sublists of arbitrary positive integers, as we're using a generalized 
    // version of the FRI protocol which can handle arbitrary folding factors
    let l_list: Vec<usize> = vec![vec![2; 11 + i]].concat();

    // We define `s` as the generator of the group of order 2^5 * n, where n is the length of the isogeny walk
    let mut s = F::new(MontFp!("20166910023067061949242007329499498203359431897882042367010355835922513064073298943410517857055490103593376746276366289375459766625"), MontFp!("9070588010778744328358501144613351095932673883221130019470787206592208103281447479286045237519441445473142536447684307414402867833"));
    let mut s = F::new(MontFp!("20166910023067061949242007329499498203359431897882042367010355835922513064073298943410517857055490103593376746276366289375459766625"), MontFp!("9070588010778744328358501144613351095932673883221130019470787206592208103281447479286045237519441445473142536447684307414402867833"));

    // We modify s based on the length of the walk
    for _ in 0..5 - i {
        s = s.pow(&[2]);
    }
    let mut g: F = s.clone();
    /// <g> is the interpolation domain for the witness polynomials
    for _ in 0..5 {
        g = g.pow(&[2]);
    }
    // Define the field element `r` as some element not in <s> such that r<s> is a coset
    let r: F = F::new(Fp::from(5), Fp::from(3));
    
    // Read the witness polynomial from a file
    let witness: DensePolynomial<F> =
        DensePolynomial { coeffs: lines_from_file(&format!("Phis/polynomial_{}.txt", 9 + i)).unwrap() };
    let n = witness.coeffs.len();

    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);

    // Create a blinding factor polynomial
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }.naive_mul(&DensePolynomial {
        coeffs: vec![vec![-F::from(1)], vec![F::from(0); n - 1], vec![F::from(1)]].concat(),
    });
    
    // Apply the blinding factor to the witness polynomial
    let b_witness: DensePolynomial<F> = witness.clone() + blinding_factor;

    // Read the psi polynomial from a file
    let psi: DensePolynomial<F> =
        DensePolynomial { coeffs: lines_from_file(&format!("Psis/polynomial_{}.txt", 9 + i)).unwrap() };

    // Evaluate the blinded witness polynomial at specific points
    let y_start: F = b_witness.evaluate(&F::from(1));
    let y_end: F = b_witness.evaluate(&g.pow(&[n as u64 - 1]));

    /// s_ord is the multiplicative order of s
    let s_ord: u64 = n as u64 * 32;
    /// The below two variables are security parameteres are for the FRI protocol
    let rep_param: usize = 1;
    let grinding_param: u8 = 32;

    /// Prove the validity of the isogeny walk
    let (challenge_vals, roots_fri, roots, paths_fri, points_fri, additional_paths_and_points, ws, paths_and_points) =
        prove(witness, psi, g, s, r, s_ord, &y_start, &y_end, l_list.clone(), rep_param, grinding_param);
    println!("Prover Time: {} s", now.elapsed().as_secs());
    let prover_time = now.elapsed().as_secs();

    let now = Instant::now();
    /// Verify the proof
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
