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
    use crate::Num;

    pub fn g_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        inlines::g_sep(s, ln_s, ln_s_pl_1)
    }

    pub fn g<T: Num>(s: T) -> T {
        inlines::g(s)
    }

    pub fn g_xi_sep<T: Num>(ln_s: T) -> T {
        inlines::g_xi_sep(ln_s)
    }

    pub fn g_xi<T: Num>(s: T) -> T {
        inlines::g_xi(s)
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
    pub extern "C" fn oneloop__ghost__g_sep(s: R, ln_s: R, ln_s_pl_1: R) -> R {
        inlines::g_sep(s, ln_s, ln_s_pl_1)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g(s: R) -> R {
        inlines::g(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g_sep__complex(s: C, ln_s: C, ln_s_pl_1: C) -> C {
        inlines::g_sep(s, ln_s, ln_s_pl_1)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g__complex(s: C) -> C {
        inlines::g(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g_xi_sep(ln_s: R) -> R {
        inlines::g_xi_sep(ln_s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g_xi(s: R) -> R {
        inlines::g_xi(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g_xi_sep__complex(ln_s: C) -> C {
        inlines::g_xi_sep(ln_s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__ghost__g_xi__complex(s: C) -> C {
        inlines::g_xi(s)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::Num;

    #[inline(always)]
    pub fn l_g<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        s_pl_1_2 * (sinv * 2. - sinv2) * ln_s_pl_1 - s * 2. * ln_s
    }

    #[inline(always)]
    pub fn r_g<T: Num>(s: T) -> T {
        s.inv() + 2.
    }

    #[inline(always)]
    pub fn g_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        (l_g(s, ln_s, ln_s_pl_1) + r_g(s)) / 12.
    }

    #[inline(always)]
    pub fn g<T: Num>(s: T) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        g_sep(s, ln_s, ln_s_pl_1)
    }

    #[inline(always)]
    pub fn g_xi_sep<T: Num>(ln_s: T) -> T {
        -ln_s / 12.
    }

    #[inline(always)]
    pub fn g_xi<T: Num>(s: T) -> T {
        -s.ln() / 12.
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::ghost;
    use crate::{Num, C, I, R};

    const TOLERANCE: R = 1e-12;

    const REAL_TEST_VAL: [R; 4] = [1., 1.834, 2.5, 8.18];
    const COMPLEX_TEST_VAL: [C; 5] = [
        C { re: 1., im: -0.5 },
        C { re: 1., im: 0.5 },
        C { re: 1.347, im: 0. },
        C {
            re: 2.56,
            im: 0.732,
        },
        C {
            re: 7.28,
            im: 5.166,
        },
    ];

    fn assert_equal<T: Num>(lhs: T, rhs: T) {
        assert!(
            (lhs - rhs).abs() < TOLERANCE,
            "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
        );
    }

    #[test]
    fn test_g() {
        use ghost::ffi::{oneloop__ghost__g, oneloop__ghost__g__complex};
        use ghost::g;

        const REAL_RESULTS: [R; 4] = [
            0.48104906018664845,
            0.579741749632399,
            0.6366840011360758,
            0.8855899339647263,
        ];
        let complex_results: [C; 5] = [
            0.4923979959662195 - 0.07124251080687459 * I,
            0.4923979959662195 + 0.07124251080687459 * I,
            0.5272889838654602 + 0. * I,
            0.6472387703649405 + 0.053746283212455476 * I,
            0.9017966251236884 + 0.14090558554076238 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(g(s), REAL_RESULTS[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__ghost__g(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(g(s), complex_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__ghost__g__complex(s), complex_results[i]));
    }

    #[test]
    fn test_g_xi() {
        use ghost::ffi::{oneloop__ghost__g_xi, oneloop__ghost__g_xi__complex};
        use ghost::g_xi;

        const REAL_RESULTS: [R; 4] = [
            -0.0,
            -0.05054161448618943,
            -0.07635756098951292,
            -0.1751410125512213,
        ];
        let complex_results: [C; 5] = [
            -0.009297647971425406 + 0.038637300750067174 * I,
            -0.009297647971425406 - 0.038637300750067174 * I,
            -0.024823324785685573 + 0. * I,
            -0.08160850400245938 - 0.023208834563292327 * I,
            -0.18242055664494833 - 0.05142917783424054 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(g_xi(s), REAL_RESULTS[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__ghost__g_xi(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(g_xi(s), complex_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__ghost__g_xi__complex(s), complex_results[i]));
    }
}
