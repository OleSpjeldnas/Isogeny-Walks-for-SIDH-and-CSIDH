use ark_ff::{Fp2, Fp2Config, Fp448, MontBackend, MontFp};

/// Define the field for the SIKE p434 parameters
#[derive(ark_ff::fp::MontConfig)]
#[modulus = "24439423661345221551909145011457493619085780243761596511325807336205221239331976725970216671828618445898719026692884939342314733567"]
#[generator = "5"]

pub struct FqConfig;
pub type F = Fp448<MontBackend<FqConfig, 7>>;

pub type Fq2 = Fp2<Fq2Config>;
pub struct Fq2Config;
pub struct Ft(Fq2);

impl Fp2Config for Fq2Config {
    type Fp = F;

    /// NONRESIDUE = -1
    const NONRESIDUE: F = MontFp!("-1");

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP2_C1: &'static [F] = &[MontFp!("-1"), MontFp!("-1")];
}
