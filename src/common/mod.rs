/// Common thermal functions.
pub mod thermal;

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

    /// Relativistic energy squared.
    ///
    /// As a function of the momentum `p` and mass `m`,
    /// `p^2 + m^2`.
    pub fn energy_squared<T: Num>(p: R, m: T) -> T {
        inlines::energy_squared(p, m)
    }

    /// Relativistic energy.
    ///
    /// As a function of the momentum `p` and mass `m`,
    /// `sqrt(p^2 + m^2)`.
    pub fn energy<T: Num>(p: R, m: T) -> T {
        inlines::energy(p, m)
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

    pub use super::thermal::ffi::*;

    #[no_mangle]
    pub extern "C" fn energy_squared(p: R, m: R) -> R {
        inlines::energy_squared(p, m)
    }

    #[no_mangle]
    pub extern "C" fn energy(p: R, m: R) -> R {
        inlines::energy(p, m)
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
    pub fn energy_squared<T: Num>(p: R, m: T) -> T {
        m * m + p * p
    }

    #[inline(always)]
    pub fn energy<T: Num>(p: R, m: T) -> T {
        energy_squared(p, m).sqrt()
    }
}
