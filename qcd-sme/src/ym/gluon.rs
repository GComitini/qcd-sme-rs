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

    pub fn dressing_inv_landau<T: Num>(s: T, f0: R) -> T {
        inlines::dressing_inv_landau(s, f0)
    }

    pub fn dressing_inv<T: Num>(s: T, f0: R, xi: R) -> T {
        inlines::dressing_inv(s, f0, xi)
    }

    pub fn dressing_landau<T: Num>(s: T, f0: R) -> T {
        inlines::dressing_landau(s, f0)
    }

    pub fn dressing<T: Num>(s: T, f0: R, xi: R) -> T {
        inlines::dressing(s, f0, xi)
    }

    pub fn propagator_landau<T: Num>(s: T, f0: R) -> T {
        inlines::propagator_landau(s, f0)
    }

    #[inline(always)]
    pub fn propagator<T: Num>(s: T, f0: R, xi: R) -> T {
        inlines::propagator(s, f0, xi)
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
    pub extern "C" fn ym__gluon__dressing_inv_landau(s: R, f0: R) -> R {
        inlines::dressing_inv_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing_inv_landau__complex(s: C, f0: R) -> C {
        inlines::dressing_inv_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing_inv(s: R, f0: R, xi: R) -> R {
        inlines::dressing_inv(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing_inv__complex(s: C, f0: R, xi: R) -> C {
        inlines::dressing_inv(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing_landau(s: R, f0: R) -> R {
        inlines::dressing_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing_landau__complex(s: C, f0: R) -> C {
        inlines::dressing_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing(s: R, f0: R, xi: R) -> R {
        inlines::dressing(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__dressing__complex(s: C, f0: R, xi: R) -> C {
        inlines::dressing(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__propagator_landau(s: R, f0: R) -> R {
        inlines::propagator_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__propagator_landau__complex(s: C, f0: R) -> C {
        inlines::propagator_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__propagator(s: R, f0: R, xi: R) -> R {
        inlines::propagator(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__gluon__propagator__complex(s: C, f0: R, xi: R) -> C {
        inlines::propagator(s, f0, xi)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::low_level::oneloop::gluon::inlines::*;
    use crate::{Num, R};

    #[inline(always)]
    pub fn dressing_inv_landau<T: Num>(s: T, f0: R) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0)
    }

    #[inline(always)]
    pub fn dressing_inv<T: Num>(s: T, f0: R, xi: R) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        ym_dressing_inv_sep(s, ln_s, ln_s_pl_1, f0, xi)
    }

    #[inline(always)]
    pub fn dressing_landau<T: Num>(s: T, f0: R) -> T {
        dressing_inv_landau(s, f0).inv()
    }

    #[inline(always)]
    pub fn dressing<T: Num>(s: T, f0: R, xi: R) -> T {
        dressing_inv(s, f0, xi).inv()
    }

    #[inline(always)]
    pub fn propagator_landau<T: Num>(s: T, f0: R) -> T {
        dressing_landau(s, f0) * s.inv()
    }

    #[inline(always)]
    pub fn propagator<T: Num>(s: T, f0: R, xi: R) -> T {
        dressing(s, f0, xi) * s.inv()
    }
}

#[cfg(test)]
mod tests {
    use crate::ym::gluon;
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
        if rhs == T::zero() {
            assert!(
                (lhs - rhs).abs() < TOLERANCE,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
            );
        } else {
            assert!(
                (lhs / rhs - 1.).abs() < TOLERANCE,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e} rhs"
            );
        }
    }

    #[test]
    fn test_dressing_inv_landau_0() {
        use gluon::dressing_inv_landau;
        use gluon::ffi::{ym__gluon__dressing_inv_landau, ym__gluon__dressing_inv_landau__complex};

        const REAL_RESULTS: [R; 4] = [
            1.25258737422535,
            1.1796181212809973,
            1.224310148974833,
            1.7199708898738115,
        ];
        let complex_results: [C; 5] = [
            1.141216547642006 + 0.10632676182751036 * I,
            1.141216547642006 - 0.10632676182751036 * I,
            1.186533662224201 + 0. * I,
            1.2225243812835327 + 0.06907180360840477 * I,
            1.740527303983705 + 0.36449585222116365 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, -0.876), REAL_RESULTS[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(ym__gluon__dressing_inv_landau(s, -0.876), REAL_RESULTS[i]);
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__gluon__dressing_inv_landau__complex(s, -0.876),
                complex_results[i],
            );
        });
    }

    fn test_dressing_inv(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use gluon::dressing_inv;
        use gluon::ffi::{ym__gluon__dressing_inv, ym__gluon__dressing_inv__complex};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, -0.876, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(ym__gluon__dressing_inv(s, -0.876, xi), real_results[i]);
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__gluon__dressing_inv__complex(s, -0.876, xi),
                complex_results[i],
            );
        });
    }

    #[test]
    fn test_dressing_inv_landau() {
        const REAL_RESULTS: [R; 4] = [
            1.25258737422535,
            1.1796181212809973,
            1.224310148974833,
            1.7199708898738115,
        ];
        let complex_results: [C; 5] = [
            1.141216547642006 + 0.10632676182751036 * I,
            1.141216547642006 - 0.10632676182751036 * I,
            1.186533662224201 + 0. * I,
            1.2225243812835327 + 0.06907180360840477 * I,
            1.740527303983705 + 0.36449585222116365 * I,
        ];
        test_dressing_inv(0., REAL_RESULTS, complex_results);
    }

    #[test]
    fn test_dressing_inv_feynman() {
        const REAL_RESULTS: [R; 4] = [
            1.3359207075586832,
            1.0886314098431225,
            1.0590022119885176,
            1.3114932289936028,
        ];
        let complex_results: [C; 5] = [
            1.171013812762538 + 0.24803303648915967 * I,
            1.171013812762538 - 0.24803303648915967 * I,
            1.1781493439750925 + 0. * I,
            1.0407408674458878 + 0.006474361373422277 * I,
            1.3124101738042617 + 0.24944668950975402 * I,
        ];
        test_dressing_inv(1., REAL_RESULTS, complex_results);
    }

    #[test]
    fn test_dressing_landau_0() {
        use gluon::dressing_landau;
        use gluon::ffi::{ym__gluon__dressing_landau, ym__gluon__dressing_landau__complex};

        let mut real_results: [R; 4] = [
            1.25258737422535,
            1.1796181212809973,
            1.224310148974833,
            1.7199708898738115,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.141216547642006 + 0.10632676182751036 * I,
            1.141216547642006 - 0.10632676182751036 * I,
            1.186533662224201 + 0. * I,
            1.2225243812835327 + 0.06907180360840477 * I,
            1.740527303983705 + 0.36449585222116365 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, -0.876), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(ym__gluon__dressing_landau(s, -0.876), real_results[i]);
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__gluon__dressing_landau__complex(s, -0.876),
                complex_results[i],
            );
        });
    }

    fn test_dressing(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use gluon::dressing;
        use gluon::ffi::{ym__gluon__dressing, ym__gluon__dressing__complex};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, -0.876, xi), real_results[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(ym__gluon__dressing(s, -0.876, xi), real_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__gluon__dressing__complex(s, -0.876, xi),
                complex_results[i],
            );
        });
    }

    #[test]
    fn test_dressing_landau() {
        let mut real_results = [
            1.25258737422535,
            1.1796181212809973,
            1.224310148974833,
            1.7199708898738115,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.141216547642006 + 0.10632676182751036 * I,
            1.141216547642006 - 0.10632676182751036 * I,
            1.186533662224201 + 0. * I,
            1.2225243812835327 + 0.06907180360840477 * I,
            1.740527303983705 + 0.36449585222116365 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }
        test_dressing(0., real_results, complex_results);
    }

    #[test]
    fn test_dressing_feynman() {
        let mut real_results = [
            1.3359207075586832,
            1.0886314098431225,
            1.0590022119885176,
            1.3114932289936028,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.171013812762538 + 0.24803303648915967 * I,
            1.171013812762538 - 0.24803303648915967 * I,
            1.1781493439750925 + 0. * I,
            1.0407408674458878 + 0.006474361373422277 * I,
            1.3124101738042617 + 0.24944668950975402 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }
        test_dressing(1., real_results, complex_results);
    }
}
