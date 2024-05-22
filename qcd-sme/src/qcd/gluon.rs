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

    pub fn dressing_inv_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        inlines::dressing_inv_landau(s, mq, f0)
    }

    pub fn dressing_inv<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        inlines::dressing_inv(s, mq, f0, xi)
    }

    pub fn dressing_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        inlines::dressing_landau(s, mq, f0)
    }

    pub fn dressing<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        inlines::dressing(s, mq, f0, xi)
    }

    pub fn dressing_crossed_inv_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        inlines::dressing_crossed_inv_landau(s, mq, f0)
    }

    pub fn dressing_crossed_inv<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        inlines::dressing_crossed_inv(s, mq, f0, xi)
    }

    pub fn dressing_crossed_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        inlines::dressing_crossed_landau(s, mq, f0)
    }

    pub fn dressing_crossed<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        inlines::dressing_crossed(s, mq, f0, xi)
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
    pub extern "C" fn qcd__gluon__dressing_inv_landau(s: R, mq: R, f0: R) -> R {
        inlines::dressing_inv_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_inv_landau__complex(s: C, mq: R, f0: R) -> C {
        inlines::dressing_inv_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_inv(s: R, mq: R, f0: R, xi: R) -> R {
        inlines::dressing_inv(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_inv__complex(s: C, mq: R, f0: R, xi: R) -> C {
        inlines::dressing_inv(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_landau(s: R, mq: R, f0: R) -> R {
        inlines::dressing_landau(s, mq, f0)
    }
    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_landau__complex(s: C, mq: R, f0: R) -> C {
        inlines::dressing_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing(s: R, mq: R, f0: R, xi: R) -> R {
        inlines::dressing(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing__complex(s: C, mq: R, f0: R, xi: R) -> C {
        inlines::dressing(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed_inv_landau(s: R, mq: R, f0: R) -> R {
        inlines::dressing_crossed_inv_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed_inv_landau__complex(s: C, mq: R, f0: R) -> C {
        inlines::dressing_crossed_inv_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed_inv(s: R, mq: R, f0: R, xi: R) -> R {
        inlines::dressing_crossed_inv(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed_inv__complex(s: C, mq: R, f0: R, xi: R) -> C {
        inlines::dressing_crossed_inv(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed_landau(s: R, mq: R, f0: R) -> R {
        inlines::dressing_crossed_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed_landau__complex(s: C, mq: R, f0: R) -> C {
        inlines::dressing_crossed_landau(s, mq, f0)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed(s: R, mq: R, f0: R, xi: R) -> R {
        inlines::dressing_crossed(s, mq, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn qcd__gluon__dressing_crossed__complex(s: C, mq: R, f0: R, xi: R) -> C {
        inlines::dressing_crossed(s, mq, f0, xi)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::consts::nf_div_nc;
    use crate::low_level::oneloop::gluon::inlines::*;
    use crate::{Num, R};

    #[inline(always)]
    pub fn dressing_inv_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        let mq2 = mq * mq;
        let s_q = s / mq2;
        ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0) + f_q(s_q) * nf_div_nc()
    }

    #[inline(always)]
    pub fn dressing_inv<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        let mq2 = mq * mq;
        let s_q = s / mq2;
        ym_dressing_inv_sep(s, ln_s, ln_s_pl_1, f0, xi) + f_q(s_q) * nf_div_nc()
    }

    #[inline(always)]
    pub fn dressing_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        dressing_inv_landau(s, mq, f0).inv()
    }

    #[inline(always)]
    pub fn dressing<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        dressing_inv(s, mq, f0, xi).inv()
    }

    #[inline(always)]
    pub fn dressing_crossed_inv_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        let mq2 = mq * mq;
        let s_q = s / mq2;
        ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0) + f_q_crossed(s_q) * nf_div_nc()
    }

    #[inline(always)]
    pub fn dressing_crossed_inv<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        let mq2 = mq * mq;
        let s_q = s / mq2;
        ym_dressing_inv_sep(s, ln_s, ln_s_pl_1, f0, xi) + f_q_crossed(s_q) * nf_div_nc()
    }

    #[inline(always)]
    pub fn dressing_crossed_landau<T: Num>(s: T, mq: R, f0: R) -> T {
        dressing_crossed_inv_landau(s, mq, f0).inv()
    }

    #[inline(always)]
    pub fn dressing_crossed<T: Num>(s: T, mq: R, f0: R, xi: R) -> T {
        dressing_crossed_inv(s, mq, f0, xi).inv()
    }
}

#[cfg(test)]
mod tests {
    use crate::qcd::gluon;
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
        if rhs != T::zero() {
            assert!(
                (lhs / rhs - 1.).abs() < TOLERANCE,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e} rhs"
            );
        } else {
            assert!(
                (lhs - rhs).abs() < TOLERANCE,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
            );
        }
    }

    #[test]
    fn test_dressing_inv_landau_0() {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_inv_landau;
        use gluon::ffi::{
            qcd__gluon__dressing_inv_landau, qcd__gluon__dressing_inv_landau__complex,
        };

        let mq = get_default_quark_mass();

        const REAL_RESULTS: [R; 4] = [
            1.2096404250687727,
            1.1138500463610683,
            1.1443964215038607,
            1.5730856514168061,
        ];
        let complex_results: [C; 5] = [
            1.096503343759402 + 0.12216298895152423 * I,
            1.096503343759402 - 0.12216298895152423 * I,
            1.1332001262908316 + 0. * I,
            1.1401568799573452 + 0.055334806887922264 * I,
            1.589666710187457 + 0.3247512836933115 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, mq, -0.876), REAL_RESULTS[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_inv_landau(s, mq, -0.876),
                REAL_RESULTS[i],
            )
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_inv_landau(s, mq, -0.876), complex_results[i])
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_inv_landau__complex(s, mq, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing_inv(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_inv;
        use gluon::ffi::{qcd__gluon__dressing_inv, qcd__gluon__dressing_inv__complex};

        let mq = get_default_quark_mass();

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, mq, -0.876, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(qcd__gluon__dressing_inv(s, mq, -0.876, xi), real_results[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, mq, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_inv__complex(s, mq, -0.876, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_inv_landau() {
        const REAL_RESULTS: [R; 4] = [
            1.2096404250687727,
            1.1138500463610683,
            1.1443964215038607,
            1.5730856514168061,
        ];
        let complex_results: [C; 5] = [
            1.096503343759402 + 0.12216298895152423 * I,
            1.096503343759402 - 0.12216298895152423 * I,
            1.1332001262908316 + 0. * I,
            1.1401568799573452 + 0.055334806887922264 * I,
            1.589666710187457 + 0.3247512836933115 * I,
        ];
        test_dressing_inv(0., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_inv_feynman() {
        const REAL_RESULTS: [R; 4] = [
            1.2929737584021062,
            1.0228633349231935,
            0.9790884845175452,
            1.1646079905365974,
        ];
        let complex_results: [C; 5] = [
            1.126300608879934 + 0.26386926361317353 * I,
            1.126300608879934 - 0.26386926361317353 * I,
            1.1248158080417232 + 0. * I,
            0.9583733661197 - 0.007262635347060228 * I,
            1.161549580008014 + 0.20970212098190188 * I,
        ];
        test_dressing_inv(1., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_landau_0() {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_landau;
        use gluon::ffi::{qcd__gluon__dressing_landau, qcd__gluon__dressing_landau__complex};

        let mq = get_default_quark_mass();

        let mut real_results: [R; 4] = [
            1.2096404250687727,
            1.1138500463610683,
            1.1443964215038607,
            1.5730856514168061,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.096503343759402 + 0.12216298895152423 * I,
            1.096503343759402 - 0.12216298895152423 * I,
            1.1332001262908316 + 0. * I,
            1.1401568799573452 + 0.055334806887922264 * I,
            1.589666710187457 + 0.3247512836933115 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, mq, -0.876), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(qcd__gluon__dressing_landau(s, mq, -0.876), real_results[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, mq, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_landau__complex(s, mq, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing;
        use gluon::ffi::{qcd__gluon__dressing, qcd__gluon__dressing__complex};

        let mq = get_default_quark_mass();

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, mq, -0.876, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(qcd__gluon__dressing(s, mq, -0.876, xi), real_results[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, mq, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing__complex(s, mq, -0.876, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_landau() {
        let mut real_results: [R; 4] = [
            1.2096404250687727,
            1.1138500463610683,
            1.1443964215038607,
            1.5730856514168061,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.096503343759402 + 0.12216298895152423 * I,
            1.096503343759402 - 0.12216298895152423 * I,
            1.1332001262908316 + 0. * I,
            1.1401568799573452 + 0.055334806887922264 * I,
            1.589666710187457 + 0.3247512836933115 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }
        test_dressing(0., real_results, complex_results)
    }

    #[test]
    fn test_dressing_feynman() {
        let mut real_results = [
            1.2929737584021062,
            1.0228633349231935,
            0.9790884845175452,
            1.1646079905365974,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.126300608879934 + 0.26386926361317353 * I,
            1.126300608879934 - 0.26386926361317353 * I,
            1.1248158080417232 + 0. * I,
            0.9583733661197 - 0.007262635347060228 * I,
            1.161549580008014 + 0.20970212098190188 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }
        test_dressing(1., real_results, complex_results)
    }

    #[test]
    fn test_dressing_crossed_inv_landau_0() {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_crossed_inv_landau;
        use gluon::ffi::{
            qcd__gluon__dressing_crossed_inv_landau,
            qcd__gluon__dressing_crossed_inv_landau__complex,
        };

        let mq = get_default_quark_mass();

        const REAL_RESULTS: [R; 4] = [
            1.1452200013339067,
            1.0278169334356733,
            1.0479069719664291,
            1.446544935725117,
        ];
        let complex_results: [C; 5] = [
            1.0283574166816276 + 0.13902433016869392 * I,
            1.0283574166816276 - 0.13902433016869392 * I,
            1.058123207942039 + 0. * I,
            1.0412407008118778 + 0.0464390633718008 * I,
            1.4593613027930297 + 0.31457224197479217 * I,
        ];

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_crossed_inv_landau(s, mq, -0.876), REAL_RESULTS[i])
        });

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed_inv_landau(s, mq, -0.876),
                REAL_RESULTS[i],
            )
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                dressing_crossed_inv_landau(s, mq, -0.876),
                complex_results[i],
            )
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed_inv_landau__complex(s, mq, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing_crossed_inv(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_crossed_inv;
        use gluon::ffi::{
            qcd__gluon__dressing_crossed_inv, qcd__gluon__dressing_crossed_inv__complex,
        };

        let mq = get_default_quark_mass();

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_crossed_inv(s, mq, -0.876, xi), real_results[i])
        });

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed_inv(s, mq, -0.876, xi),
                real_results[i],
            )
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_crossed_inv(s, mq, -0.876, xi), complex_results[i])
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed_inv__complex(s, mq, -0.876, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_crossed_inv_landau() {
        const REAL_RESULTS: [R; 4] = [
            1.1452200013339067,
            1.0278169334356733,
            1.0479069719664291,
            1.446544935725117,
        ];
        let complex_results: [C; 5] = [
            1.0283574166816276 + 0.13902433016869392 * I,
            1.0283574166816276 - 0.13902433016869392 * I,
            1.058123207942039 + 0. * I,
            1.0412407008118778 + 0.0464390633718008 * I,
            1.4593613027930297 + 0.31457224197479217 * I,
        ];
        test_dressing_crossed_inv(0., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_crossed_inv_feynman() {
        const REAL_RESULTS: [R; 4] = [
            1.22855333466724,
            0.9368302219977984,
            0.8825990349801137,
            1.0380672748449082,
        ];
        let complex_results: [C; 5] = [
            1.0581546818021597 + 0.28073060483034323 * I,
            1.0581546818021597 - 0.28073060483034323 * I,
            1.0497388896929305 + 0. * I,
            0.8594571869742331 - 0.016158378863181694 * I,
            1.0312441726135864 + 0.19952307926338256 * I,
        ];
        test_dressing_crossed_inv(1., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_crossed_landau_0() {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_crossed_landau;
        use gluon::ffi::{
            qcd__gluon__dressing_crossed_landau, qcd__gluon__dressing_crossed_landau__complex,
        };

        let mq = get_default_quark_mass();

        let mut real_results: [R; 4] = [
            1.1452200013339067,
            1.0278169334356733,
            1.0479069719664291,
            1.446544935725117,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.0283574166816276 + 0.13902433016869392 * I,
            1.0283574166816276 - 0.13902433016869392 * I,
            1.058123207942039 + 0. * I,
            1.0412407008118778 + 0.0464390633718008 * I,
            1.4593613027930297 + 0.31457224197479217 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_crossed_landau(s, mq, -0.876), real_results[i])
        });

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed_landau(s, mq, -0.876),
                real_results[i],
            )
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_crossed_landau(s, mq, -0.876), complex_results[i])
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed_landau__complex(s, mq, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing_crossed(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use crate::consts::get_default_quark_mass;
        use gluon::dressing_crossed;
        use gluon::ffi::{qcd__gluon__dressing_crossed, qcd__gluon__dressing_crossed__complex};

        let mq = get_default_quark_mass();

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_crossed(s, mq, -0.876, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed(s, mq, -0.876, xi),
                real_results[i],
            )
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(dressing_crossed(s, mq, -0.876, xi), complex_results[i])
        });

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                qcd__gluon__dressing_crossed__complex(s, mq, -0.876, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_crossed_landau() {
        let mut real_results: [R; 4] = [
            1.1452200013339067,
            1.0278169334356733,
            1.0479069719664291,
            1.446544935725117,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.0283574166816276 + 0.13902433016869392 * I,
            1.0283574166816276 - 0.13902433016869392 * I,
            1.058123207942039 + 0. * I,
            1.0412407008118778 + 0.0464390633718008 * I,
            1.4593613027930297 + 0.31457224197479217 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }
        test_dressing_crossed(0., real_results, complex_results)
    }

    #[test]
    fn test_dressing_crossed_feynman() {
        let mut real_results = [
            1.22855333466724,
            0.9368302219977984,
            0.8825990349801137,
            1.0380672748449082,
        ];
        for v in &mut real_results {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.0581546818021597 + 0.28073060483034323 * I,
            1.0581546818021597 - 0.28073060483034323 * I,
            1.0497388896929305 + 0. * I,
            0.8594571869742331 - 0.016158378863181694 * I,
            1.0312441726135864 + 0.19952307926338256 * I,
        ];
        for v in &mut complex_results {
            *v = v.inv();
        }
        test_dressing_crossed(1., real_results, complex_results)
    }
}
