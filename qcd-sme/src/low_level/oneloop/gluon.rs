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

    pub fn f_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        inlines::f_sep(s, ln_s, ln_s_pl_1)
    }

    pub fn f<T: Num>(s: T) -> T {
        inlines::f(s)
    }

    pub fn f_xi_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        inlines::f_xi_sep(s, ln_s, ln_s_pl_1)
    }

    pub fn f_xi<T: Num>(s: T) -> T {
        inlines::f_xi(s)
    }

    pub fn f_q<T: Num>(s: T) -> T {
        inlines::f_q(s)
    }

    pub fn f_q_crossed<T: Num>(s: T) -> T {
        inlines::f_q_crossed(s)
    }

    pub fn ym_dressing_inv_landau_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T, f0: R) -> T {
        inlines::ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0)
    }

    pub fn ym_dressing_inv_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T, f0: R, xi: R) -> T {
        inlines::ym_dressing_inv_sep(s, ln_s, ln_s_pl_1, f0, xi)
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
    pub extern "C" fn oneloop__gluon__f_sep(s: R, ln_s: R, ln_s_pl_1: R) -> R {
        inlines::f_sep(s, ln_s, ln_s_pl_1)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f(s: R) -> R {
        inlines::f(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_sep__complex(s: C, ln_s: C, ln_s_pl_1: C) -> C {
        inlines::f_sep(s, ln_s, ln_s_pl_1)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f__complex(s: C) -> C {
        inlines::f(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi_sep(s: R, ln_s: R, ln_s_pl_1: R) -> R {
        inlines::f_xi_sep(s, ln_s, ln_s_pl_1)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi(s: R) -> R {
        inlines::f_xi(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi_sep__complex(s: C, ln_s: C, ln_s_pl_1: C) -> C {
        inlines::f_xi_sep(s, ln_s, ln_s_pl_1)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi__complex(s: C) -> C {
        inlines::f_xi(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_q(s: R) -> R {
        inlines::f_q(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_q__complex(s: C) -> C {
        inlines::f_q(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_q_crossed(s: R) -> R {
        inlines::f_q_crossed(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_q_crossed__complex(s: C) -> C {
        inlines::f_q_crossed(s)
    }

    #[no_mangle]
    pub extern "C" fn ym__oneloop__gluon__dressing_inv_landau_sep(
        s: R,
        ln_s: R,
        ln_s_pl_1: R,
        f0: R,
    ) -> R {
        inlines::ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__oneloop__gluon__dressing_inv_landau_sep__complex(
        s: C,
        ln_s: C,
        ln_s_pl_1: C,
        f0: R,
    ) -> C {
        inlines::ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0)
    }

    #[no_mangle]
    pub extern "C" fn ym__oneloop__gluon__dressing_inv_sep(
        s: R,
        ln_s: R,
        ln_s_pl_1: R,
        f0: R,
        xi: R,
    ) -> R {
        inlines::ym_dressing_inv_sep(s, ln_s, ln_s_pl_1, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ym__oneloop__gluon__dressing_inv_sep__complex(
        s: C,
        ln_s: C,
        ln_s_pl_1: C,
        f0: R,
        xi: R,
    ) -> C {
        inlines::ym_dressing_inv_sep(s, ln_s, ln_s_pl_1, f0, xi)
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
    pub fn l_a<T: Num>(s: T, ln_s: T) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let a = (sinv * 4. + 1.).sqrt();
        (s2 * 3. - s * 34. - 28. - sinv * 24.) * a * (((a - 1.) / 2.).ln() * 2. + ln_s)
    }

    #[inline(always)]
    pub fn l_b<T: Num>(s: T, ln_s_pl_1: T) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let sinv3 = sinv * sinv2;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        s_pl_1_2 * (-sinv * 20. + 3. + sinv2 * 11. - sinv3 * 2.) * 2. * ln_s_pl_1
    }

    #[inline(always)]
    pub fn l_c<T: Num>(s: T, ln_s: T) -> T {
        let s2 = s * s;
        (-s2 * 3. + 2.) * ln_s
    }

    #[inline(always)]
    pub fn r_a<T: Num>(s: T) -> T {
        let sinv = s.inv();
        -(s + 4.) * (s - 20. + sinv * 12.)
    }

    #[inline(always)]
    pub fn r_b<T: Num>(s: T) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        s_pl_1_2 * 2. * (-sinv * 10. + 1. + sinv2)
    }

    #[inline(always)]
    pub fn r_c<T: Num>(s: T) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        sinv2 * 2. + 2. - s2
    }

    #[inline(always)]
    pub fn f_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        let sinv = s.inv();
        (l_a(s, ln_s) + l_b(s, ln_s_pl_1) + l_c(s, ln_s) + r_a(s) + r_b(s) + r_c(s)) / 72.
            + sinv * 5. / 8.
    }

    #[inline(always)]
    pub fn f<T: Num>(s: T) -> T {
        let s_pl_1 = s + 1.;
        let ln_s = s.ln();
        let ln_s_pl_1 = s_pl_1.ln();
        f_sep(s, ln_s, ln_s_pl_1)
    }

    #[inline(always)]
    pub fn f_xi_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        sinv / 4.
            - (s * 2. * ln_s - (sinv2 - sinv) * (sinv - s2) * 2. * ln_s_pl_1 + 3. - sinv * 3.
                + sinv2 * 2.)
                / 12.
    }

    #[inline(always)]
    pub fn f_xi<T: Num>(s: T) -> T {
        let ln_s = s.ln();
        let s_pl_1 = s + 1.;
        let ln_s_pl_1 = s_pl_1.ln();
        f_xi_sep(s, ln_s, ln_s_pl_1)
    }

    // Defined so that f_q(0)=0.
    //
    // Don't define a "sep" variant as s is p^2/M_QUARK^2 and its logarithms
    // are unlikely to appear in other expressions.
    #[inline(always)]
    pub fn f_q<T: Num>(s: T) -> T {
        let sinv = s.inv();
        let a2 = sinv * 4. + 1.;
        let a = a2.sqrt();
        ((-sinv + 0.5) * a * ((a - 1.) / (a + 1.)).ln() - sinv * 2. + 5. / 6.) * 4. / 9.
    }

    // Equal to f_q(s)+2s*df_q(s)/ds modulo a constant defined so that f_q_crossed(0)=0.
    //
    // Don't define a "sep" variant as s is p^2/M_QUARK^2 and its logarithms
    // are unlikely to appear in other expressions.
    #[inline(always)]
    pub fn f_q_crossed<T: Num>(s: T) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let a2 = sinv * 4. + 1.;
        let a = a2.sqrt();
        (sinv * 4. - 1. / 6. + (sinv + 0.5 + sinv2 * 8.) / a * ((a - 1.) / (a + 1.)).ln()) * 4. / 9.
    }

    #[inline(always)]
    pub fn ym_dressing_inv_landau_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T, f0: R) -> T {
        f_sep(s, ln_s, ln_s_pl_1) + f0
    }

    #[inline(always)]
    pub fn ym_dressing_inv_sep<T: Num>(s: T, ln_s: T, ln_s_pl_1: T, f0: R, xi: R) -> T {
        ym_dressing_inv_landau_sep(s, ln_s, ln_s_pl_1, f0) + f_xi_sep(s, ln_s, ln_s_pl_1) * xi
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::gluon;
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
    fn test_f() {
        use gluon::f;
        use gluon::ffi::{oneloop__gluon__f, oneloop__gluon__f__complex};

        const REAL_RESULTS: [R; 4] = [
            2.12858737422535,
            2.055618121280997,
            2.100310148974833,
            2.5959708898738114,
        ];
        let complex_results: [C; 5] = [
            2.017216547642006 + 0.10632676182751036 * I,
            2.017216547642006 - 0.10632676182751036 * I,
            2.062533662224201 + 0. * I,
            2.0985243812835326 + 0.06907180360840477 * I,
            2.616527303983705 + 0.36449585222116365 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f(s), REAL_RESULTS[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f(s), complex_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f__complex(s), complex_results[i]));
    }

    #[test]
    fn test_f_xi() {
        use gluon::f_xi;
        use gluon::ffi::{oneloop__gluon__f_xi, oneloop__gluon__f_xi__complex};

        const REAL_RESULTS: [R; 4] = [
            0.08333333333333334,
            -0.09098671143787487,
            -0.1653079369863152,
            -0.40847766088020865,
        ];
        let complex_results: [C; 5] = [
            0.029797265120532046 + 0.1417062746616493 * I,
            0.029797265120532046 - 0.1417062746616493 * I,
            -0.008384318249108513 + 0. * I,
            -0.181783513837645 - 0.0625974422349825 * I,
            -0.42811713017944325 - 0.11504916271140962 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_xi(s), REAL_RESULTS[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f_xi(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_xi(s), complex_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f_xi__complex(s), complex_results[i]));
    }

    #[test]
    fn test_f_q() {
        use gluon::f_q;
        use gluon::ffi::{oneloop__gluon__f_q, oneloop__gluon__f_q__complex};

        const REAL_RESULTS: [R; 4] = [
            -0.040286361891847346,
            -0.06882604993530607,
            -0.0891399017500785,
            -0.2116655881071791,
        ];
        let complex_results: [C; 5] = [
            -0.04108195144887958 + 0.018275406354830614 * I,
            -0.04108195144887958 - 0.018275406354830614 * I,
            -0.05263306134948051 + 0. * I,
            -0.09192015626362243 - 0.021049245646939316 * I,
            -0.21395381256611867 - 0.08636454168837476 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_q(s), REAL_RESULTS[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f_q(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_q(s), complex_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f_q__complex(s), complex_results[i]));
    }

    #[test]
    fn test_f_q_crossed() {
        use gluon::f_q_crossed;
        use gluon::ffi::{oneloop__gluon__f_q_crossed, oneloop__gluon__f_q_crossed__complex};

        const REAL_RESULTS: [R; 4] = [
            -0.11357849147764694,
            -0.1862365230440458,
            -0.23458327877629956,
            -0.48104183141541057,
        ];
        let complex_results: [C; 5] = [
            -0.11697385059874044 + 0.04829797592211286 * I,
            -0.11697385059874044 - 0.04829797592211286 * I,
            -0.14573365482653244 + 0. * I,
            -0.24231647920553134 - 0.04836593654467078 * I,
            -0.49722997194803575 - 0.14978025966031802 * I,
        ];

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_q_crossed(s), REAL_RESULTS[i]));

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(oneloop__gluon__f_q_crossed(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_q_crossed(s), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(oneloop__gluon__f_q_crossed__complex(s), complex_results[i])
        });
    }
}
