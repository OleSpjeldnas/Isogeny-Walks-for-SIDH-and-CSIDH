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

/// Implement FftField for Fq2
/// This is used for FFT-based interpolation and multiplication
use ark_ff::FftField;

impl FftField for Fq2 {
    // The generator of the multiplicative group of the field
    const GENERATOR: Self = Fq2::new(MontFp!("2"), MontFp!("1"));


    // The largest power of two dividing (q - 1).
    const TWO_ADICITY: u32 = 217;

    // A 2^TWO_ADICITY-th primitive root of unity
    // For example, some element g where g^(2^TWO_ADICITY) = 1,
    // but g^(2^(TWO_ADICITY-1)) != 1, etc.
    const TWO_ADIC_ROOT_OF_UNITY: Self = Fq2::new(MontFp!("1265641318332283385680345622106262255016986152776970677011496306402637098857925419143401023483759358905674114613032730660970188006"), MontFp!("24302508882629253666484860493871125283132764152316444477996794014562063591984863827165794089737987196964883650291668738548192491206"));
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{FftField, Field}; // Bring in the traits

    // We'll need BigInteger for exponentiation
    use ark_ff::BigInteger;
    use ark_std::test_rng; // optional if you want to do random tests

    #[test]
    fn test_generator_is_not_one() {
        // Just a sanity check that Fq2::GENERATOR != 1
        let g = Fq2::GENERATOR;
        assert_ne!(g, Fq2::ONE, "Generator must not be 1");
    }

    #[test]
    fn test_two_adic_root_of_unity_order() {
        let root = Fq2::TWO_ADIC_ROOT_OF_UNITY;

        // 1) Check root^(2^TWO_ADICITY) == 1
        // Build the exponent 2^217 as a BigInteger.
        // For TWO_ADICITY = 217, that's 1 << 217.
        // We'll use `<<` on a `u64` then convert to the field's BigInteger type.

        let exp_2_to_k = {
            let mut e = <F::BigInt as From<u64>>::from(1u64); 
            e = e << Fq2::TWO_ADICITY; // shift left by 217 bits
            e
        };

        let val = root.pow(exp_2_to_k);
        assert_eq!(val, Fq2::ONE, "root^(2^TWO_ADICITY) must be 1");

        // 2) (Optional) Check root^(2^(k - 1)) != 1 to confirm it's *primitive*
        // i.e., it has *exactly* order 2^k, not a divisor.
        if Fq2::TWO_ADICITY > 0 {
            let exp_2_to_k_minus_1 = {
                let mut e = <F::BigInt as From<u64>>::from(1u64);
                e = e << (Fq2::TWO_ADICITY - 1);
                e
            };
            let val_smaller = root.pow(exp_2_to_k_minus_1);
            assert_ne!(val_smaller, Fq2::ONE, 
                "root should have order 2^k, so raising it to 2^(k-1) should not give 1"
            );
        }
    }
}

