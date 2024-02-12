use num::complex::ComplexFloat;

/// A real number.
pub type R = f64;
/// A complex number.
pub type C = num::Complex<R>;
/// An integer type for the number of colors.
///
/// Used by [set_number_of_colors](crate::consts::set_number_of_colors).
pub type NCTYPE = u32;
/// An integer type for the number of fermions.
///
/// Used by [set_number_of_fermions](crate::consts::set_number_of_fermions).
pub type NFTYPE = u32;

/// A number.
pub trait Num:
    Copy
    + num::Num
    + std::fmt::Debug
    + std::fmt::Display
    + std::ops::Add<R, Output = Self>
    + std::ops::Mul<R, Output = Self>
    + std::ops::Div<R, Output = Self>
    + std::ops::Sub<R, Output = Self>
    + std::ops::Neg<Output = Self>
{
    fn abs(&self) -> R;
    fn exp(&self) -> Self;
    fn inv(&self) -> Self;
    fn ln(&self) -> Self;
    fn sqrt(&self) -> Self;
}

impl Num for R {
    #[inline(always)]
    fn abs(&self) -> R {
        (*self as R).abs()
    }

    #[inline(always)]
    fn exp(&self) -> Self {
        (*self as R).exp()
    }

    #[inline(always)]
    fn inv(&self) -> Self {
        1. / self
    }

    #[inline(always)]
    fn ln(&self) -> Self {
        (*self as R).ln()
    }

    #[inline(always)]
    fn sqrt(&self) -> Self {
        (*self as R).sqrt()
    }
}

impl Num for C {
    #[inline(always)]
    fn abs(&self) -> R {
        (*self as C).abs()
    }

    #[inline(always)]
    fn exp(&self) -> Self {
        (*self as C).exp()
    }

    #[inline(always)]
    fn inv(&self) -> Self {
        1. / self
    }

    #[inline(always)]
    fn ln(&self) -> Self {
        (*self as C).ln()
    }

    #[inline(always)]
    fn sqrt(&self) -> Self {
        (*self as C).sqrt()
    }
}
