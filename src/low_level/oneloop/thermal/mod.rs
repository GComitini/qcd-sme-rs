pub mod ghost;
pub mod gluon;

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
    use crate::{Num, C, R};

    pub fn r_0<T1: Num, T2: Num>(om: T1, p: R, en: T2, delta: R) -> C {
        inlines::r_0(om, p, en, delta)
    }

    pub fn tlog<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2, delta: R) -> C {
        inlines::tlog(q, om, p, en, delta)
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
    use crate::{C, R};

    #[no_mangle]
    pub extern "C" fn oneloop__tlog(q: R, om: R, p: R, en: R, delta: R) -> C {
        inlines::tlog(q, om, p, en, delta)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::{Num, C, I, R};

    #[inline(always)]
    pub fn r_0<T1: Num, T2: Num>(om: T1, p: R, en: T2, delta: R) -> C {
        om * om + om * (en * 2. * I) + p * p + delta
    }

    #[inline(always)]
    pub fn tlog<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2, delta: R) -> C {
        ((r_0(om, p, en, delta) + q * 2. * p) / (r_0(om, p, en, delta) - q * 2. * p)).ln()
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::thermal;
    use crate::{Num, C, I, R};

    const TOLERANCE: R = 1e-12;

    fn assert_equal<T: Num>(lhs: T, rhs: T) {
        assert!(
            (lhs - rhs).abs() < TOLERANCE,
            "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
        );
    }

    #[test]
    fn test_tlog() {
        use thermal::{ffi::oneloop__tlog, tlog};

        let args: [(R, R, R, R, R); 6] = [
            (0.03, 0., 2.12, 4., 2.),
            (0.03, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 1.27, 4., 2.),
            (0.047, 0.18, 1.27, 6.9, 2.),
            (0.047, 0.18, 1.27, 6.9, 3.21),
        ];

        let res: [C; 6] = [
            0.039177220079546576 + 0. * I,
            0.03717215443439166 - 0.008203228369481758 * I,
            0.05824470539037811 - 0.012858080275916684 * I,
            0.05666490383705108 - 0.022398146165630385 * I,
            0.0447257555087438 - 0.03049215271192196 * I,
            0.03897509874778243 - 0.01994626073948082 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(tlog(*q, *om, *p, *en, *delta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(oneloop__tlog(*q, *om, *p, *en, *delta), res[i])
            });
    }
}
