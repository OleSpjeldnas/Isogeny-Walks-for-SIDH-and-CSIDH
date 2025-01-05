use ark_ff::{Fp2, Fp2Config, Fp448, MontBackend, MontFp, One, UniformRand, Zero};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize, CanonicalSerializeWithFlags, Flags, Read,
    SerializationError, Write,
};
use rand::Rng;
use std::fmt::{self, Debug, Display};
use std::hash::Hash;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use zeroize::Zeroize;
/// Define the field for the SIKE p434 parameters
#[derive(ark_ff::fp::MontConfig)]
#[modulus = "24439423661345221551909145011457493619085780243761596511325807336205221239331976725970216671828618445898719026692884939342314733567"]
#[generator = "5"]

pub struct FqConfig;
pub type F = Fp448<MontBackend<FqConfig, 7>>;

pub type Fq2 = Fp2<Fq2Config>;
pub struct Fq2Config;

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

impl FftField for MyFq2 {
    // The generator of the multiplicative group of the field
    const GENERATOR: Self = Self(Fq2::new(MontFp!("2"), MontFp!("1")));

    // The largest power of two dividing (q - 1).
    const TWO_ADICITY: u32 = 217;

    // A 2^TWO_ADICITY-th primitive root of unity
    // For example, some element g where g^(2^TWO_ADICITY) = 1,
    // but g^(2^(TWO_ADICITY-1)) != 1, etc.
    const TWO_ADIC_ROOT_OF_UNITY: Self = Self(Fq2::new(MontFp!("1265641318332283385680345622106262255016986152776970677011496306402637098857925419143401023483759358905674114613032730660970188006"), MontFp!("24302508882629253666484860493871125283132764152316444477996794014562063591984863827165794089737987196964883650291668738548192491206")));
}

#[derive(
    Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Zeroize, CanonicalSerialize, CanonicalDeserialize,
)]
pub struct MyFq2(pub Fq2);

// Implement Display
impl Display for MyFq2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// Implement Zero and One
impl Zero for MyFq2 {
    fn zero() -> Self {
        Self(Fq2::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for MyFq2 {
    fn one() -> Self {
        Self(Fq2::one())
    }

    fn is_one(&self) -> bool {
        self.0.is_one()
    }
}

// Implement UniformRand
impl UniformRand for MyFq2 {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self(Fq2::rand(rng))
    }
}

// Implement arithmetic operations
impl Add for MyFq2 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}

impl Sub for MyFq2 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl Neg for MyFq2 {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}

impl Mul for MyFq2 {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}

impl Div for MyFq2 {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

// Implement arithmetic assign operations
impl AddAssign for MyFq2 {
    fn add_assign(&mut self, other: Self) {
        self.0 += other.0;
    }
}

impl SubAssign for MyFq2 {
    fn sub_assign(&mut self, other: Self) {
        self.0 -= other.0;
    }
}

impl MulAssign for MyFq2 {
    fn mul_assign(&mut self, other: Self) {
        self.0 *= other.0;
    }
}

impl DivAssign for MyFq2 {
    fn div_assign(&mut self, other: Self) {
        self.0 /= other.0;
    }
}

// Implement reference operations
impl<'a> Add<&'a MyFq2> for MyFq2 {
    type Output = Self;
    fn add(self, other: &'a MyFq2) -> Self {
        Self(self.0 + other.0)
    }
}

impl<'a> Sub<&'a MyFq2> for MyFq2 {
    type Output = Self;
    fn sub(self, other: &'a MyFq2) -> Self {
        Self(self.0 - other.0)
    }
}

impl<'a> Mul<&'a MyFq2> for MyFq2 {
    type Output = Self;
    fn mul(self, other: &'a MyFq2) -> Self {
        Self(self.0 * other.0)
    }
}

impl<'a> Div<&'a MyFq2> for MyFq2 {
    type Output = Self;
    fn div(self, other: &'a MyFq2) -> Self {
        Self(self.0 / other.0)
    }
}

// Implement reference assign operations
impl<'a> AddAssign<&'a MyFq2> for MyFq2 {
    fn add_assign(&mut self, other: &'a MyFq2) {
        self.0 += other.0;
    }
}

impl<'a> SubAssign<&'a MyFq2> for MyFq2 {
    fn sub_assign(&mut self, other: &'a MyFq2) {
        self.0 -= other.0;
    }
}

impl<'a> MulAssign<&'a MyFq2> for MyFq2 {
    fn mul_assign(&mut self, other: &'a MyFq2) {
        self.0 *= other.0;
    }
}

impl<'a> DivAssign<&'a MyFq2> for MyFq2 {
    fn div_assign(&mut self, other: &'a MyFq2) {
        self.0 /= other.0;
    }
}

// Implement Sum and Product
impl Sum for MyFq2 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a> Sum<&'a MyFq2> for MyFq2 {
    fn sum<I: Iterator<Item = &'a MyFq2>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl Product for MyFq2 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a> Product<&'a MyFq2> for MyFq2 {
    fn product<I: Iterator<Item = &'a MyFq2>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

use ark_ff::{fields::Field, LegendreSymbol, SqrtPrecomputation};

impl Field for MyFq2 {
    type BasePrimeField = F;
    type BasePrimeFieldIter = std::array::IntoIter<F, 2>;

    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = None;
    const ONE: Self = Self(Fq2::ONE);
    const ZERO: Self = Self(Fq2::ZERO);

    fn extension_degree() -> u64 {
        2
    }

    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter {
        [self.0.c0, self.0.c1].into_iter()
    }

    fn from_base_prime_field_elems(elems: &[Self::BasePrimeField]) -> Option<Self> {
        if elems.len() != 2 {
            return None;
        }
        Some(Self(Fq2::new(elems[0], elems[1])))
    }

    fn from_base_prime_field(elem: Self::BasePrimeField) -> Self {
        Self(Fq2::new(elem, Self::BasePrimeField::zero()))
    }

    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        Fq2::from_random_bytes_with_flags(bytes).map(|(f, flags)| (Self(f), flags))
    }

    fn legendre(&self) -> LegendreSymbol {
        self.0.legendre()
    }

    fn square(&self) -> Self {
        Self(self.0.square())
    }

    fn square_in_place(&mut self) -> &mut Self {
        self.0.square_in_place();
        self
    }

    fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Self)
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        if let Some(inv) = self.0.inverse() {
            self.0 = inv;
            Some(self)
        } else {
            None
        }
    }

    fn double(&self) -> Self {
        Self(self.0.double())
    }

    fn double_in_place(&mut self) -> &mut Self {
        self.0.double_in_place();
        self
    }

    fn neg_in_place(&mut self) -> &mut Self {
        self.0.neg_in_place();
        self
    }

    fn frobenius_map(&mut self, power: usize) {
        self.0.frobenius_map(power);
    }
}

// Implement From<T> for various integer types
macro_rules! impl_from_int {
    ($($t:ty),*) => {
        $(
            impl From<$t> for MyFq2 {
                fn from(value: $t) -> Self {
                    Self(Fq2::from(value))
                }
            }
        )*
    };
}

impl_from_int!(u128, u64, u32, u16, u8, i128, i64, i32, i16, i8, bool);

// Implement division operations for mutable references
impl<'a> Div<&'a mut MyFq2> for MyFq2 {
    type Output = Self;
    fn div(self, other: &'a mut MyFq2) -> Self {
        Self(self.0 / other.0)
    }
}

impl<'a> DivAssign<&'a mut MyFq2> for MyFq2 {
    fn div_assign(&mut self, other: &'a mut MyFq2) {
        self.0 /= other.0;
    }
}

// Implement arithmetic operations with mutable references
impl<'a> Add<&'a mut MyFq2> for MyFq2 {
    type Output = Self;
    fn add(self, other: &'a mut MyFq2) -> Self {
        Self(self.0 + other.0)
    }
}

impl<'a> Sub<&'a mut MyFq2> for MyFq2 {
    type Output = Self;
    fn sub(self, other: &'a mut MyFq2) -> Self {
        Self(self.0 - other.0)
    }
}

impl<'a> Mul<&'a mut MyFq2> for MyFq2 {
    type Output = Self;
    fn mul(self, other: &'a mut MyFq2) -> Self {
        Self(self.0 * other.0)
    }
}

// Implement arithmetic assign operations with mutable references
impl<'a> AddAssign<&'a mut MyFq2> for MyFq2 {
    fn add_assign(&mut self, other: &'a mut MyFq2) {
        self.0 += other.0;
    }
}

impl<'a> SubAssign<&'a mut MyFq2> for MyFq2 {
    fn sub_assign(&mut self, other: &'a mut MyFq2) {
        self.0 -= other.0;
    }
}

impl<'a> MulAssign<&'a mut MyFq2> for MyFq2 {
    fn mul_assign(&mut self, other: &'a mut MyFq2) {
        self.0 *= other.0;
    }
}

// Implement CanonicalSerializeWithFlags
impl CanonicalSerializeWithFlags for MyFq2 {
    fn serialize_with_flags<W: Write, F: Flags>(&self, mut writer: W, flags: F) -> Result<(), SerializationError> {
        self.0.serialize_with_flags(&mut writer, flags)
    }

    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        self.0.serialized_size_with_flags::<F>()
    }
}

// Implement CanonicalDeserializeWithFlags
impl CanonicalDeserializeWithFlags for MyFq2 {
    fn deserialize_with_flags<R: Read, F: Flags>(mut reader: R) -> Result<(Self, F), SerializationError> {
        let (f, flags) = Fq2::deserialize_with_flags(&mut reader)?;
        Ok((MyFq2(f), flags))
    }
}

impl MyFq2 {
    /// Creates a new field element from two base field elements
    pub fn new(c0: F, c1: F) -> Self {
        Self(Fq2::new(c0, c1))
    }

    /// Returns the first coefficient (c0)
    pub fn c0(&self) -> F {
        self.0.c0
    }

    /// Returns the second coefficient (c1)
    pub fn c1(&self) -> F {
        self.0.c1
    }
}

// Implement subtraction for references
impl<'a, 'b> Sub<&'b MyFq2> for &'a MyFq2 {
    type Output = MyFq2;
    fn sub(self, other: &'b MyFq2) -> MyFq2 {
        MyFq2(self.0 - other.0)
    }
}

// Also implement multiplication for references for completeness
impl<'a, 'b> Mul<&'b MyFq2> for &'a MyFq2 {
    type Output = MyFq2;
    fn mul(self, other: &'b MyFq2) -> MyFq2 {
        MyFq2(self.0 * other.0)
    }
}

// And addition for references
impl<'a, 'b> Add<&'b MyFq2> for &'a MyFq2 {
    type Output = MyFq2;
    fn add(self, other: &'b MyFq2) -> MyFq2 {
        MyFq2(self.0 + other.0)
    }
}
