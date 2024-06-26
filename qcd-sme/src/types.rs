use num::complex::ComplexFloat;
use std::fmt;
use std::ops;

/// An integration method.
pub use peroxide::numerical::integral::Integral;

/// A real number.
pub type R = f64;

/// A complex number.
pub type C = num::Complex<R>;

/// An integer type for the number of colors.
///
/// Used by [`set_number_of_colors`](crate::consts::set_number_of_colors).
pub type NCTYPE = u32;

/// An integer type for the number of fermions.
///
/// Used by [`set_number_of_fermions`](crate::consts::set_number_of_fermions).
pub type NFTYPE = u32;

/// A number.
pub trait Num:
    Copy
    + num::Num
    + ops::AddAssign
    + ops::SubAssign
    + ops::MulAssign
    + ops::DivAssign
    + fmt::Debug
    + fmt::Display
    + ops::Add<R, Output = Self>
    + ops::Mul<R, Output = Self>
    + ops::Div<R, Output = Self>
    + ops::Sub<R, Output = Self>
    + ops::AddAssign<R>
    + ops::SubAssign<R>
    + ops::MulAssign<R>
    + ops::DivAssign<R>
    + ops::Neg<Output = Self>
    + ops::Add<C, Output = C>
    + ops::Sub<C, Output = C>
    + ops::Mul<C, Output = C>
    + ops::Div<C, Output = C>
    + From<R>
    + Into<C>
{
    fn abs(&self) -> R;
    fn exp(&self) -> Self;
    fn im(&self) -> R;
    fn inv(&self) -> Self;
    fn ln(&self) -> Self;
    fn re(&self) -> R;
    fn sqrt(&self) -> Self;
}

impl Num for R {
    #[inline(always)]
    fn abs(&self) -> R {
        (*self as Self).abs()
    }

    #[inline(always)]
    fn exp(&self) -> Self {
        (*self as Self).exp()
    }

    #[inline(always)]
    fn im(&self) -> R {
        0.
    }

    #[inline(always)]
    fn inv(&self) -> Self {
        1. / self
    }

    #[inline(always)]
    fn ln(&self) -> Self {
        (*self as Self).ln()
    }

    #[inline(always)]
    fn re(&self) -> R {
        *self
    }

    #[inline(always)]
    fn sqrt(&self) -> Self {
        (*self as Self).sqrt()
    }
}

impl Num for C {
    #[inline(always)]
    fn abs(&self) -> R {
        (*self as Self).abs()
    }

    #[inline(always)]
    fn exp(&self) -> Self {
        (*self as Self).exp()
    }

    #[inline(always)]
    fn im(&self) -> R {
        self.im
    }

    #[inline(always)]
    fn inv(&self) -> Self {
        1. / self
    }

    #[inline(always)]
    fn ln(&self) -> Self {
        (*self as Self).ln()
    }

    #[inline(always)]
    fn re(&self) -> R {
        self.re
    }

    #[inline(always)]
    fn sqrt(&self) -> Self {
        (*self as Self).sqrt()
    }
}
