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

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_inv_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
    ) -> T {
        inlines::dressing_inv_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    pub fn dressing_inv_landau<T: Num>(s: T, f0: R) -> T {
        inlines::dressing_inv_landau(s, f0)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_inv_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
        xi: R,
    ) -> T {
        inlines::dressing_inv_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    pub fn dressing_inv<T: Num>(s: T, f0: R, xi: R) -> T {
        inlines::dressing_inv(s, f0, xi)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
    ) -> T {
        inlines::dressing_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    pub fn dressing_landau<T: Num>(s: T, f0: R) -> T {
        inlines::dressing_landau(s, f0)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
        xi: R,
    ) -> T {
        inlines::dressing_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    pub fn dressing<T: Num>(s: T, f0: R, xi: R) -> T {
        inlines::dressing(s, f0, xi)
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
    pub extern "C" fn gluon__dressing_inv_landau_sep(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        s_q: R,
        f0: R,
    ) -> R {
        inlines::dressing_inv_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv_landau(s: R, f0: R) -> R {
        inlines::dressing_inv_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv_landau_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        s_q: C,
        f0: R,
    ) -> C {
        inlines::dressing_inv_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv_landau__complex(s: C, f0: R) -> C {
        inlines::dressing_inv_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv_sep(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        s_q: R,
        f0: R,
        xi: R,
    ) -> R {
        inlines::dressing_inv_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv(s: R, f0: R, xi: R) -> R {
        inlines::dressing_inv(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        s_q: C,
        f0: R,
        xi: R,
    ) -> C {
        inlines::dressing_inv_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_inv__complex(s: C, f0: R, xi: R) -> C {
        inlines::dressing_inv(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_landau_sep(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        s_q: R,
        f0: R,
    ) -> R {
        inlines::dressing_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_landau(s: R, f0: R) -> R {
        inlines::dressing_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_landau_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        s_q: C,
        f0: R,
    ) -> C {
        inlines::dressing_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_landau__complex(s: C, f0: R) -> C {
        inlines::dressing_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_sep(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        s_q: R,
        f0: R,
        xi: R,
    ) -> R {
        inlines::dressing_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing(s: R, f0: R, xi: R) -> R {
        inlines::dressing(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        s_q: C,
        f0: R,
        xi: R,
    ) -> C {
        inlines::dressing_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn gluon__dressing__complex(s: C, f0: R, xi: R) -> C {
        inlines::dressing(s, f0, xi)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::consts::{m_quark, nf_div_nc};
    use crate::low_level::oneloop::gluon::inlines::*;
    use crate::{Num, R};

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_inv_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
    ) -> T {
        f_sep(s2, s, sinv, sinv2, s_pl_1_2) + f_q(s_q) * nf_div_nc() + f0
    }

    #[inline(always)]
    pub fn dressing_inv_landau<T: Num>(s: T, f0: R) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        let x = m_quark();
        let x2 = x * x;
        let s_q = s / x2;
        dressing_inv_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_inv_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
        xi: R,
    ) -> T {
        dressing_inv_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0)
            + f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2) * xi
    }

    #[inline(always)]
    pub fn dressing_inv<T: Num>(s: T, f0: R, xi: R) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        let x = m_quark();
        let x2 = x * x;
        let s_q = s / x2;
        dressing_inv_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi)
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
    ) -> T {
        dressing_inv_landau_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0).inv()
    }

    #[inline(always)]
    pub fn dressing_landau<T: Num>(s: T, f0: R) -> T {
        dressing_inv_landau(s, f0).inv()
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        s_q: T,
        f0: R,
        xi: R,
    ) -> T {
        dressing_inv_sep(s2, s, sinv, sinv2, s_pl_1_2, s_q, f0, xi).inv()
    }

    #[inline(always)]
    pub fn dressing<T: Num>(s: T, f0: R, xi: R) -> T {
        dressing_inv(s, f0, xi).inv()
    }
}

#[cfg(test)]
mod tests {
    use crate::gluon;
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
    fn test_dressing_inv_landau_0() {
        use gluon::dressing_inv_landau;
        use gluon::ffi::{gluon__dressing_inv_landau, gluon__dressing_inv_landau__complex};

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

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, -0.876), REAL_RESULTS[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(gluon__dressing_inv_landau(s, -0.876), REAL_RESULTS[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                gluon__dressing_inv_landau__complex(s, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing_inv(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use gluon::dressing_inv;
        use gluon::ffi::{gluon__dressing_inv, gluon__dressing_inv__complex};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, -0.876, xi), real_results[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(gluon__dressing_inv(s, -0.876, xi), real_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                gluon__dressing_inv__complex(s, -0.876, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_inv_landau() {
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
        test_dressing_inv(0., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_inv_feynman() {
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
        test_dressing_inv(1., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_landau_0() {
        use gluon::dressing_landau;
        use gluon::ffi::{gluon__dressing_landau, gluon__dressing_landau__complex};

        let mut real_results: [R; 4] = [
            1.1452200013339067,
            1.0278169334356733,
            1.0479069719664291,
            1.446544935725117,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.0283574166816276 + 0.13902433016869392 * I,
            1.0283574166816276 - 0.13902433016869392 * I,
            1.058123207942039 + 0. * I,
            1.0412407008118778 + 0.0464390633718008 * I,
            1.4593613027930297 + 0.31457224197479217 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, -0.876), real_results[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(gluon__dressing_landau(s, -0.876), real_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                gluon__dressing_landau__complex(s, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use gluon::dressing;
        use gluon::ffi::{gluon__dressing, gluon__dressing__complex};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, -0.876, xi), real_results[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(gluon__dressing(s, -0.876, xi), real_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(gluon__dressing__complex(s, -0.876, xi), complex_results[i])
        });
    }

    #[test]
    fn test_dressing_landau() {
        let mut real_results: [R; 4] = [
            1.1452200013339067,
            1.0278169334356733,
            1.0479069719664291,
            1.446544935725117,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.0283574166816276 + 0.13902433016869392 * I,
            1.0283574166816276 - 0.13902433016869392 * I,
            1.058123207942039 + 0. * I,
            1.0412407008118778 + 0.0464390633718008 * I,
            1.4593613027930297 + 0.31457224197479217 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }
        test_dressing(0., real_results, complex_results)
    }

    #[test]
    fn test_dressing_feynman() {
        let mut real_results = [
            1.22855333466724,
            0.9368302219977984,
            0.8825990349801137,
            1.0380672748449082,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.0581546818021597 + 0.28073060483034323 * I,
            1.0581546818021597 - 0.28073060483034323 * I,
            1.0497388896929305 + 0. * I,
            0.8594571869742331 - 0.016158378863181694 * I,
            1.0312441726135864 + 0.19952307926338256 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }
        test_dressing(1., real_results, complex_results)
    }
}