pub use native::*;

// For use in Rust
//
// - Collect into a module to improve code organization, but immediately re-export.
//
// - Each function simply calls its inline analogue with the sole objective of not
//   inlining code when not necessary. The inline functions are still available to
//   the crate for more specialized use.
//
// - Make pub(crate) only in this module so it can be re-exported by crate::ghost.
pub(crate) mod native {
    use super::inlines;
    use crate::{Num, R};

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_inv_landau_sep<T: Num>(
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        ln_s: T,
        g0: R,
    ) -> T {
        inlines::dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    pub fn dressing_inv_landau<T: Num>(s: T, g0: R) -> T {
        inlines::dressing_inv_landau(s, g0)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_inv_sep<T: Num>(
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        ln_s: T,
        g0: R,
        xi: R,
    ) -> T {
        inlines::dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    pub fn dressing_inv<T: Num>(s: T, g0: R, xi: R) -> T {
        inlines::dressing_inv(s, g0, xi)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_landau_sep<T: Num>(s: T, sinv: T, sinv2: T, s_pl_1_2: T, ln_s: T, g0: R) -> T {
        inlines::dressing_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    pub fn dressing_landau<T: Num>(s: T, g0: R) -> T {
        inlines::dressing_landau(s, g0)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_sep<T: Num>(s: T, sinv: T, sinv2: T, s_pl_1_2: T, ln_s: T, g0: R, xi: R) -> T {
        inlines::dressing_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    pub fn dressing<T: Num>(s: T, g0: R, xi: R) -> T {
        inlines::dressing(s, g0, xi)
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
    pub extern "C" fn ym__ghost__dressing_inv_landau_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        g0: R,
    ) -> R {
        inlines::dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv_landau(s: R, g0: R) -> R {
        inlines::dressing_inv_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv_landau_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        g0: R,
    ) -> C {
        inlines::dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv_landau__complex(s: C, g0: R) -> C {
        inlines::dressing_inv_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        g0: R,
        xi: R,
    ) -> R {
        inlines::dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv(s: R, g0: R, xi: R) -> R {
        inlines::dressing_inv(s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        g0: R,
        xi: R,
    ) -> C {
        inlines::dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_inv__complex(s: C, g0: R, xi: R) -> C {
        inlines::dressing_inv(s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_landau_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        g0: R,
    ) -> R {
        inlines::dressing_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_landau(s: R, g0: R) -> R {
        inlines::dressing_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_landau_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        g0: R,
    ) -> C {
        inlines::dressing_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_landau__complex(s: C, g0: R) -> C {
        inlines::dressing_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        g0: R,
        xi: R,
    ) -> R {
        inlines::dressing_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing(s: R, g0: R, xi: R) -> R {
        inlines::dressing(s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        g0: R,
        xi: R,
    ) -> C {
        inlines::dressing_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__ghost__dressing__complex(s: C, g0: R, xi: R) -> C {
        inlines::dressing(s, g0, xi)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::low_level::oneloop::ghost::inlines::*;
    use crate::{Num, R};

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_inv_landau_sep<T: Num>(
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        ln_s: T,
        g0: R,
    ) -> T {
        g_sep(s, sinv, sinv2, s_pl_1_2, ln_s) + g0
    }

    #[inline(always)]
    pub fn dressing_inv_landau<T: Num>(s: T, g0: R) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        let ln_s = s.ln();
        dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0)
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_inv_sep<T: Num>(
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_1_2: T,
        ln_s: T,
        g0: R,
        xi: R,
    ) -> T {
        dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0) + g_xi_sep(ln_s) * xi
    }

    #[inline(always)]
    pub fn dressing_inv<T: Num>(s: T, g0: R, xi: R) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        let ln_s = s.ln();
        dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi)
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_landau_sep<T: Num>(s: T, sinv: T, sinv2: T, s_pl_1_2: T, ln_s: T, g0: R) -> T {
        dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0).inv()
    }

    #[inline(always)]
    pub fn dressing_landau<T: Num>(s: T, g0: R) -> T {
        dressing_inv_landau(s, g0).inv()
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_sep<T: Num>(s: T, sinv: T, sinv2: T, s_pl_1_2: T, ln_s: T, g0: R, xi: R) -> T {
        dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, g0, xi).inv()
    }

    #[inline(always)]
    pub fn dressing<T: Num>(s: T, g0: R, xi: R) -> T {
        dressing_inv(s, g0, xi).inv()
    }
}

#[cfg(test)]
mod tests {
    use crate::ym::ghost;
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
        use ghost::dressing_inv_landau;
        use ghost::ffi::{ym__ghost__dressing_inv_landau, ym__ghost__dressing_inv_landau__complex};

        const REAL_RESULTS: [R; 4] = [
            0.6262890601866484,
            0.724981749632399,
            0.7819240011360759,
            1.0308299339647262,
        ];
        let complex_results: [C; 5] = [
            0.6376379959662195 - 0.07124251080687459 * I,
            0.6376379959662195 + 0.07124251080687459 * I,
            0.6725289838654602 + 0. * I,
            0.7924787703649405 + 0.053746283212455476 * I,
            1.0470366251236884 + 0.14090558554076238 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, 0.14524), REAL_RESULTS[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(ym__ghost__dressing_inv_landau(s, 0.14524), REAL_RESULTS[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, 0.14524), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__ghost__dressing_inv_landau__complex(s, 0.14524),
                complex_results[i],
            )
        });
    }

    fn test_dressing_inv(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use ghost::dressing_inv;
        use ghost::ffi::{ym__ghost__dressing_inv, ym__ghost__dressing_inv__complex};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, 0.14524, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(ym__ghost__dressing_inv(s, 0.14524, xi), real_results[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, 0.14524, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__ghost__dressing_inv__complex(s, 0.14524, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_inv_landau() {
        const REAL_RESULTS: [R; 4] = [
            0.6262890601866484,
            0.724981749632399,
            0.7819240011360759,
            1.0308299339647262,
        ];
        let complex_results: [C; 5] = [
            0.6376379959662195 - 0.07124251080687459 * I,
            0.6376379959662195 + 0.07124251080687459 * I,
            0.6725289838654602 + 0. * I,
            0.7924787703649405 + 0.053746283212455476 * I,
            1.0470366251236884 + 0.14090558554076238 * I,
        ];
        test_dressing_inv(0., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_inv_feynman() {
        const REAL_RESULTS: [R; 4] = [
            0.6262890601866484,
            0.6744401351462096,
            0.7055664401465629,
            0.8556889214135049,
        ];
        let complex_results: [C; 5] = [
            0.6283403479947941 - 0.032605210056807415 * I,
            0.6283403479947941 + 0.032605210056807415 * I,
            0.6477056590797746 + 0. * I,
            0.7108702663624812 + 0.03053744864916315 * I,
            0.8646160684787401 + 0.08947640770652185 * I,
        ];
        test_dressing_inv(1., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_landau_0() {
        use ghost::dressing_landau;
        use ghost::ffi::{ym__ghost__dressing_landau, ym__ghost__dressing_landau__complex};

        let mut real_results: [R; 4] = [
            0.6262890601866484,
            0.724981749632399,
            0.7819240011360759,
            1.0308299339647262,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            0.6376379959662195 - 0.07124251080687459 * I,
            0.6376379959662195 + 0.07124251080687459 * I,
            0.6725289838654602 + 0. * I,
            0.7924787703649405 + 0.053746283212455476 * I,
            1.0470366251236884 + 0.14090558554076238 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, 0.14524), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(ym__ghost__dressing_landau(s, 0.14524), real_results[i])
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, 0.14524), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__ghost__dressing_landau__complex(s, 0.14524),
                complex_results[i],
            )
        });
    }

    fn test_dressing(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use ghost::dressing;
        use ghost::ffi::{ym__ghost__dressing, ym__ghost__dressing__complex};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, 0.14524, xi), real_results[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(ym__ghost__dressing(s, 0.14524, xi), real_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, 0.14524, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                ym__ghost__dressing__complex(s, 0.14524, xi),
                complex_results[i],
            )
        });
    }

    #[test]
    fn test_dressing_landau() {
        let mut real_results = [
            0.6262890601866484,
            0.724981749632399,
            0.7819240011360759,
            1.0308299339647262,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            0.6376379959662195 - 0.07124251080687459 * I,
            0.6376379959662195 + 0.07124251080687459 * I,
            0.6725289838654602 + 0. * I,
            0.7924787703649405 + 0.053746283212455476 * I,
            1.0470366251236884 + 0.14090558554076238 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }
        test_dressing(0., real_results, complex_results)
    }

    #[test]
    fn test_dressing_feynman() {
        let mut real_results = [
            0.6262890601866484,
            0.6744401351462096,
            0.7055664401465629,
            0.8556889214135049,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            0.6283403479947941 - 0.032605210056807415 * I,
            0.6283403479947941 + 0.032605210056807415 * I,
            0.6477056590797746 + 0. * I,
            0.7108702663624812 + 0.03053744864916315 * I,
            0.8646160684787401 + 0.08947640770652185 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }
        test_dressing(1., real_results, complex_results)
    }
}
