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

    #[allow(clippy::too_many_arguments)]
    pub fn f_sep<T: Num>(s2: T, s: T, sinv: T, sinv2: T, s_pl_1_2: T) -> T {
        inlines::f_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    pub fn f<T: Num>(s: T) -> T {
        inlines::f(s)
    }

    pub fn f_xi_sep<T: Num>(s2: T, s: T, sinv: T, sinv2: T, s_pl_1_2: T) -> T {
        inlines::f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    pub fn f_xi<T: Num>(s: T) -> T {
        inlines::f_xi(s)
    }

    pub fn f_q<T: Num>(s: T) -> T {
        inlines::f_q(s)
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
    pub extern "C" fn oneloop__gluon__f_sep(s2: R, s: R, sinv: R, sinv2: R, s_pl_1_2: R) -> R {
        inlines::f_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f(s: R) -> R {
        inlines::f(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
    ) -> C {
        inlines::f_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f__complex(s: C) -> C {
        inlines::f(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi_sep(s2: R, s: R, sinv: R, sinv2: R, s_pl_1_2: R) -> R {
        inlines::f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi(s: R) -> R {
        inlines::f_xi(s)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__f_xi_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
    ) -> C {
        inlines::f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2)
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
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::Num;

    #[inline(always)]
    pub fn l_a<T: Num>(s2: T, s: T, sinv: T, ln_s: T) -> T {
        let a = (sinv * 4. + 1.).sqrt();
        (s2 * 3. - s * 34. - 28. - sinv * 24.) * a * (((a - 1.) / 2.).ln() * 2. + ln_s)
    }

    #[inline(always)]
    pub fn l_b<T: Num>(sinv: T, sinv2: T, s_pl_1_2: T) -> T {
        let sinv3 = sinv * sinv2;
        s_pl_1_2 * (-sinv * 20. + 3. + sinv2 * 11. - sinv3 * 2.) * s_pl_1_2.ln()
    }

    #[inline(always)]
    pub fn l_c<T: Num>(s2: T, ln_s: T) -> T {
        (-s2 * 3. + 2.) * ln_s
    }

    #[inline(always)]
    pub fn r_a<T: Num>(s: T, sinv: T) -> T {
        -(s + 4.) * (s - 20. + sinv * 12.)
    }

    #[inline(always)]
    pub fn r_b<T: Num>(sinv: T, sinv2: T, s_pl_1_2: T) -> T {
        s_pl_1_2 * 2. * (-sinv * 10. + 1. + sinv2)
    }

    #[inline(always)]
    pub fn r_c<T: Num>(s2: T, sinv2: T) -> T {
        sinv2 * 2. + 2. - s2
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn f_sep<T: Num>(s2: T, s: T, sinv: T, sinv2: T, s_pl_1_2: T) -> T {
        let ln_s = s.ln();
        (l_a(s2, s, sinv, ln_s)
            + l_b(sinv, sinv2, s_pl_1_2)
            + l_c(s2, ln_s)
            + r_a(s, sinv)
            + r_b(sinv, sinv2, s_pl_1_2)
            + r_c(s2, sinv2))
            / 72.
            + sinv * 5. / 8.
    }

    #[inline(always)]
    pub fn f<T: Num>(s: T) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        f_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    #[inline(always)]
    pub fn f_xi_sep<T: Num>(s2: T, s: T, sinv: T, sinv2: T, s_pl_1_2: T) -> T {
        sinv / 4.
            - (s * 2. * s.ln() - (sinv2 - sinv) * (sinv - s2) * s_pl_1_2.ln() + 3. - sinv * 3.
                + sinv2 * 2.)
                / 12.
    }

    #[inline(always)]
    pub fn f_xi<T: Num>(s: T) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    // Defined so that f_q(0)=0.
    //
    // Don't define a "sep" variant as s is p^2/M_QUARK^2 and its powers
    // are unlikely to appear in other expressions.
    #[inline(always)]
    pub fn f_q<T: Num>(s: T) -> T {
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let a2 = sinv * 4. + 1.;
        let a = a2.sqrt();
        (sinv * 4. - 1. / 6. + (sinv + 0.5 + sinv2 * 8.) / a * ((a - 1.) / (a + 1.)).ln()) * 4. / 9.
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
        assert!(
            (lhs - rhs).abs() < TOLERANCE,
            "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
        );
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
}