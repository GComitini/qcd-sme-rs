pub use native::*;

// For use in Rust
//
// - Collect into a module to improve code organization, but immediately re-export.
//
// - Each function simply calls its inline analogue with the sole objective of not
//   inlining code when not necessary. The inline functions are still available to
//   the crate for more specialized use.
mod native {
    use super::inlines;
    use crate::{Num, R};

    pub fn fermi_distribution<T1: Num, T2: Num>(en: T1, beta: R, mu: T2) -> T1
    where
        T1: std::ops::Sub<T2, Output = T1>,
    {
        inlines::fermi_distribution(en, beta, mu)
    }

    pub fn fermi_distribution_double<T1: Num, T2: Num>(en: T1, beta: R, mu: T2) -> T1
    where
        T1: std::ops::Sub<T2, Output = T1>,
    {
        inlines::fermi_distribution_double(en, beta, mu)
    }

    pub fn fermi_distribution_zero_chempot<T: Num>(en: T, beta: R) -> T {
        inlines::fermi_distribution_zero_chempot(en, beta)
    }

    pub fn bose_distribution<T1: Num, T2: Num>(en: T1, beta: R, mu: T2) -> T1
    where
        T1: std::ops::Sub<T2, Output = T1>,
    {
        inlines::bose_distribution(en, beta, mu)
    }

    pub fn bose_distribution_zero_chempot<T: Num>(en: T, beta: R) -> T {
        inlines::bose_distribution_zero_chempot(en, beta)
    }
}

// For use in other languages, e.g. C/C++/Python
//
// - Re-export at crate::ffi, since symbols need to be unmangled anyway and the
//   namespace will not be preserved.
//
// - Here too each function calls its inline analogue, but the objective is to
//   switch from generic to concrete argument types so that the functions can
//   be compiled into a C dynamic library. To do so, we need to double their
//   number (one function for real arguments, another for complex arguments).
pub(crate) mod ffi {
    use super::inlines;
    use crate::R;

    #[no_mangle]
    pub extern "C" fn fermi_distribution(en: R, beta: R, mu: R) -> R {
        inlines::fermi_distribution(en, beta, mu)
    }

    #[no_mangle]
    pub extern "C" fn fermi_distribution_double(en: R, beta: R, mu: R) -> R {
        inlines::fermi_distribution_double(en, beta, mu)
    }

    #[no_mangle]
    pub extern "C" fn fermi_distribution_zero_chempot(en: R, beta: R) -> R {
        inlines::fermi_distribution_zero_chempot(en, beta)
    }

    #[no_mangle]
    pub extern "C" fn bose_distribution(en: R, beta: R, mu: R) -> R {
        inlines::bose_distribution(en, beta, mu)
    }

    #[no_mangle]
    pub extern "C" fn bose_distribution_zero_chempot(en: R, beta: R) -> R {
        inlines::bose_distribution_zero_chempot(en, beta)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub mod inlines {
    use crate::low_level::common::inlines::{
        statistic_distribution_exponential, statistic_distribution_exponential_zero_chempot,
    };
    use crate::{Num, R};

    #[inline(always)]
    pub fn fermi_distribution<T1: Num, T2: Num>(en: T1, beta: R, mu: T2) -> T1
    where
        T1: std::ops::Sub<T2, Output = T1>,
    {
        (statistic_distribution_exponential(en, beta, mu) + 1.).inv()
    }

    #[inline(always)]
    pub fn fermi_distribution_double<T1: Num, T2: Num>(en: T1, beta: R, mu: T2) -> T1
    where
        T1: std::ops::Sub<T2, Output = T1>,
    {
        (statistic_distribution_exponential(en, beta, mu) + 1.).inv()
            + (statistic_distribution_exponential(en, beta, -mu) + 1.).inv()
    }

    #[inline(always)]
    pub fn fermi_distribution_zero_chempot<T: Num>(en: T, beta: R) -> T {
        (statistic_distribution_exponential_zero_chempot(en, beta) + 1.).inv()
    }

    #[inline(always)]
    pub fn bose_distribution<T1: Num, T2: Num>(en: T1, beta: R, mu: T2) -> T1
    where
        T1: std::ops::Sub<T2, Output = T1>,
    {
        (statistic_distribution_exponential(en, beta, mu) - 1.).inv()
    }

    #[inline(always)]
    pub fn bose_distribution_zero_chempot<T: Num>(en: T, beta: R) -> T {
        (statistic_distribution_exponential_zero_chempot(en, beta) - 1.).inv()
    }
}
