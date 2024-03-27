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

    pub fn statistic_distribution_exponential<T1: Num + std::ops::Sub<T2, Output = T1>, T2: Num>(
        en: T1,
        beta: R,
        mu: T2,
    ) -> T1 {
        inlines::statistic_distribution_exponential(en, beta, mu)
    }

    pub fn statistic_distribution_exponential_zero_chempot<T: Num>(en: T, beta: R) -> T {
        inlines::statistic_distribution_exponential_zero_chempot(en, beta)
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
    pub extern "C" fn statistic_distribution_exponential(en: R, beta: R, mu: R) -> R {
        inlines::statistic_distribution_exponential(en, beta, mu)
    }

    #[no_mangle]
    pub extern "C" fn statistic_distribution_exponential_zero_chempot(en: R, beta: R) -> R {
        inlines::statistic_distribution_exponential_zero_chempot(en, beta)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::{Num, R};

    #[inline(always)]
    pub fn statistic_distribution_exponential<T1: Num + std::ops::Sub<T2, Output = T1>, T2: Num>(
        en: T1,
        beta: R,
        mu: T2,
    ) -> T1 {
        ((en - mu) * beta).exp()
    }

    #[inline(always)]
    pub fn statistic_distribution_exponential_zero_chempot<T: Num>(en: T, beta: R) -> T {
        (en * beta).exp()
    }
}
