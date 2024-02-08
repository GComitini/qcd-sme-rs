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
    pub fn f_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
    ) -> T {
        inlines::f_sep(s2, s, sinv, sinv2, s_pl_4, sqrt_s_pl_4, sqrt_s, s_pl_1_2)
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

    #[allow(clippy::too_many_arguments)]
    pub fn dressing_inv_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
    ) -> T {
        inlines::dressing_inv_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
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
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
        xi: R,
    ) -> T {
        inlines::dressing_inv_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
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
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
    ) -> T {
        inlines::dressing_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
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
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
        xi: R,
    ) -> T {
        inlines::dressing_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
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
    pub extern "C" fn one_loop__gluon__f_sep__real(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_4: R,
        sqrt_s_pl_4: R,
        sqrt_s: R,
        s_pl_1_2: R,
    ) -> R {
        inlines::f_sep(s2, s, sinv, sinv2, s_pl_4, sqrt_s_pl_4, sqrt_s, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f__real(s: R) -> R {
        inlines::f(s)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_4: C,
        sqrt_s_pl_4: C,
        sqrt_s: C,
        s_pl_1_2: C,
    ) -> C {
        inlines::f_sep(s2, s, sinv, sinv2, s_pl_4, sqrt_s_pl_4, sqrt_s, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f__complex(s: C) -> C {
        inlines::f(s)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f_xi_sep__real(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
    ) -> R {
        inlines::f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f_xi__real(s: R) -> R {
        inlines::f_xi(s)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f_xi_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
    ) -> C {
        inlines::f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__f_xi__complex(s: C) -> C {
        inlines::f_xi(s)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv_landau_sep__real(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_4: R,
        sqrt_s_pl_4: R,
        sqrt_s: R,
        s_pl_1_2: R,
        f0: R,
    ) -> R {
        inlines::dressing_inv_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv_landau__real(s: R, f0: R) -> R {
        inlines::dressing_inv_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv_landau_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_4: C,
        sqrt_s_pl_4: C,
        sqrt_s: C,
        s_pl_1_2: C,
        f0: R,
    ) -> C {
        inlines::dressing_inv_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv_landau__complex(s: C, f0: R) -> C {
        inlines::dressing_inv_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv_sep__real(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_4: R,
        sqrt_s_pl_4: R,
        sqrt_s: R,
        s_pl_1_2: R,
        f0: R,
        xi: R,
    ) -> R {
        inlines::dressing_inv_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv__real(s: R, f0: R, xi: R) -> R {
        inlines::dressing_inv(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_4: C,
        sqrt_s_pl_4: C,
        sqrt_s: C,
        s_pl_1_2: C,
        f0: R,
        xi: R,
    ) -> C {
        inlines::dressing_inv_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_inv__complex(s: C, f0: R, xi: R) -> C {
        inlines::dressing_inv(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_landau_sep__real(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_4: R,
        sqrt_s_pl_4: R,
        sqrt_s: R,
        s_pl_1_2: R,
        f0: R,
    ) -> R {
        inlines::dressing_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_landau__real(s: R, f0: R) -> R {
        inlines::dressing_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_landau_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_4: C,
        sqrt_s_pl_4: C,
        sqrt_s: C,
        s_pl_1_2: C,
        f0: R,
    ) -> C {
        inlines::dressing_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_landau__complex(s: C, f0: R) -> C {
        inlines::dressing_landau(s, f0)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_sep__real(
        s2: R,
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_4: R,
        sqrt_s_pl_4: R,
        sqrt_s: R,
        s_pl_1_2: R,
        f0: R,
        xi: R,
    ) -> R {
        inlines::dressing_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing__real(s: R, f0: R, xi: R) -> R {
        inlines::dressing(s, f0, xi)
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing_sep__complex(
        s2: C,
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_4: C,
        sqrt_s_pl_4: C,
        sqrt_s: C,
        s_pl_1_2: C,
        f0: R,
        xi: R,
    ) -> C {
        inlines::dressing_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
    }

    #[no_mangle]
    pub extern "C" fn one_loop__gluon__dressing__complex(s: C, f0: R, xi: R) -> C {
        inlines::dressing(s, f0, xi)
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
    pub fn l_a<T: Num>(s2: T, s: T, sinv: T, sqrt_s_pl_4: T, sqrt_s: T) -> T {
        (s2 * 3. - s * 34. - 28. - sinv * 24.) * sqrt_s_pl_4 / sqrt_s
            * ((sqrt_s_pl_4 - sqrt_s) / (sqrt_s_pl_4 + sqrt_s)).ln()
    }

    #[inline(always)]
    pub fn l_b<T: Num>(s2: T, s: T, sinv: T, s_pl_1_2: T) -> T {
        s_pl_1_2 * (s2 * 3. - s * 20. + 11. - sinv * 2.) / s2 * s_pl_1_2.ln()
    }

    #[inline(always)]
    pub fn l_c<T: Num>(s2: T, s: T) -> T {
        (-s2 * 3. + 2.) * s.ln()
    }

    #[inline(always)]
    pub fn r_a<T: Num>(s: T, sinv: T, s_pl_4: T) -> T {
        -s_pl_4 * (s - 20. + sinv * 12.)
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
    pub fn f_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
    ) -> T {
        (l_a(s2, s, sinv, sqrt_s_pl_4, sqrt_s)
            + l_b(s2, s, sinv, s_pl_1_2)
            + l_c(s2, s)
            + r_a(s, sinv, s_pl_4)
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
        let s_pl_4 = s + 4.;
        let sqrt_s_pl_4 = s_pl_4.sqrt();
        let sqrt_s = s.sqrt();
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        f_sep(s2, s, sinv, sinv2, s_pl_4, sqrt_s_pl_4, sqrt_s, s_pl_1_2)
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

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_inv_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
    ) -> T {
        f_sep(s2, s, sinv, sinv2, s_pl_4, sqrt_s_pl_4, sqrt_s, s_pl_1_2) + f0
    }

    #[inline(always)]
    pub fn dressing_inv_landau<T: Num>(s: T, f0: R) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_4 = s + 4.;
        let sqrt_s_pl_4 = s_pl_4.sqrt();
        let sqrt_s = s.sqrt();
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        dressing_inv_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_inv_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
        xi: R,
    ) -> T {
        dressing_inv_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        ) + f_xi_sep(s2, s, sinv, sinv2, s_pl_1_2) * xi
    }

    #[inline(always)]
    pub fn dressing_inv<T: Num>(s: T, f0: R, xi: R) -> T {
        let s2 = s * s;
        let sinv = s.inv();
        let sinv2 = sinv * sinv;
        let s_pl_4 = s + 4.;
        let sqrt_s_pl_4 = s_pl_4.sqrt();
        let sqrt_s = s.sqrt();
        let s_pl_1 = s + 1.;
        let s_pl_1_2 = s_pl_1 * s_pl_1;
        dressing_inv_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub fn dressing_landau_sep<T: Num>(
        s2: T,
        s: T,
        sinv: T,
        sinv2: T,
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
    ) -> T {
        dressing_inv_landau_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
        )
        .inv()
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
        s_pl_4: T,
        sqrt_s_pl_4: T,
        sqrt_s: T,
        s_pl_1_2: T,
        f0: R,
        xi: R,
    ) -> T {
        dressing_inv_sep(
            s2,
            s,
            sinv,
            sinv2,
            s_pl_4,
            sqrt_s_pl_4,
            sqrt_s,
            s_pl_1_2,
            f0,
            xi,
        )
        .inv()
    }

    #[inline(always)]
    pub fn dressing<T: Num>(s: T, f0: R, xi: R) -> T {
        dressing_inv(s, f0, xi).inv()
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::one_loop::gluon;
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
        use gluon::ffi::{one_loop__gluon__f__complex, one_loop__gluon__f__real};

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
            .for_each(|(i, &s)| assert_equal(one_loop__gluon__f__real(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f(s), complex_results[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(one_loop__gluon__f__complex(s), complex_results[i]));
    }

    #[test]
    fn test_f_xi() {
        use gluon::f_xi;
        use gluon::ffi::{one_loop__gluon__f_xi__complex, one_loop__gluon__f_xi__real};

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
            .for_each(|(i, &s)| assert_equal(one_loop__gluon__f_xi__real(s), REAL_RESULTS[i]));

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(f_xi(s), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(one_loop__gluon__f_xi__complex(s), complex_results[i])
        });
    }

    #[test]
    fn test_dressing_inv_landau_0() {
        use gluon::dressing_inv_landau;
        use gluon::ffi::{
            one_loop__gluon__dressing_inv_landau__complex,
            one_loop__gluon__dressing_inv_landau__real,
        };

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
            assert_equal(
                one_loop__gluon__dressing_inv_landau__real(s, -0.876),
                REAL_RESULTS[i],
            )
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv_landau(s, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing_inv_landau__complex(s, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing_inv(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use gluon::dressing_inv;
        use gluon::ffi::{
            one_loop__gluon__dressing_inv__complex, one_loop__gluon__dressing_inv__real,
        };

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, -0.876, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing_inv__real(s, -0.876, xi),
                real_results[i],
            )
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_inv(s, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing_inv__complex(s, -0.876, xi),
                complex_results[i],
            )
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
        test_dressing_inv(0., REAL_RESULTS, complex_results)
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
        test_dressing_inv(1., REAL_RESULTS, complex_results)
    }

    #[test]
    fn test_dressing_landau_0() {
        use gluon::dressing_landau;
        use gluon::ffi::{
            one_loop__gluon__dressing_landau__complex, one_loop__gluon__dressing_landau__real,
        };

        let mut real_results: [R; 4] = [
            1.25258737422535,
            1.1796181212809973,
            1.224310148974833,
            1.7199708898738115,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.141216547642006 + 0.10632676182751036 * I,
            1.141216547642006 - 0.10632676182751036 * I,
            1.186533662224201 + 0. * I,
            1.2225243812835327 + 0.06907180360840477 * I,
            1.740527303983705 + 0.36449585222116365 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, -0.876), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing_landau__real(s, -0.876),
                real_results[i],
            )
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing_landau(s, -0.876), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing_landau__complex(s, -0.876),
                complex_results[i],
            )
        });
    }

    fn test_dressing(xi: R, real_results: [R; 4], complex_results: [C; 5]) {
        use gluon::dressing;
        use gluon::ffi::{one_loop__gluon__dressing__complex, one_loop__gluon__dressing__real};

        REAL_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, -0.876, xi), real_results[i]));

        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing__real(s, -0.876, xi),
                real_results[i],
            )
        });

        COMPLEX_TEST_VAL
            .iter()
            .enumerate()
            .for_each(|(i, &s)| assert_equal(dressing(s, -0.876, xi), complex_results[i]));

        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
            assert_equal(
                one_loop__gluon__dressing__complex(s, -0.876, xi),
                complex_results[i],
            )
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
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.141216547642006 + 0.10632676182751036 * I,
            1.141216547642006 - 0.10632676182751036 * I,
            1.186533662224201 + 0. * I,
            1.2225243812835327 + 0.06907180360840477 * I,
            1.740527303983705 + 0.36449585222116365 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }
        test_dressing(0., real_results, complex_results)
    }

    #[test]
    fn test_dressing_feynman() {
        let mut real_results = [
            1.3359207075586832,
            1.0886314098431225,
            1.0590022119885176,
            1.3114932289936028,
        ];
        for v in real_results.iter_mut() {
            *v = v.inv();
        }
        let mut complex_results: [C; 5] = [
            1.171013812762538 + 0.24803303648915967 * I,
            1.171013812762538 - 0.24803303648915967 * I,
            1.1781493439750925 + 0. * I,
            1.0407408674458878 + 0.006474361373422277 * I,
            1.3124101738042617 + 0.24944668950975402 * I,
        ];
        for v in complex_results.iter_mut() {
            *v = v.inv();
        }
        test_dressing(1., real_results, complex_results)
    }
}
