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
    use peroxide::numerical::integral::Integral;

    pub fn polarization_glue_l_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::polarization_glue_l_thermal_part_landau_i(q, om, p, m, beta)
    }

    pub fn polarization_glue_t_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::polarization_glue_t_thermal_part_landau_i(q, om, p, m, beta)
    }

    pub fn polarization_quark_l_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_landau_i(q, om, p, mq, beta, mu)
    }

    pub fn polarization_quark_t_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_landau_i(q, om, p, mq, beta, mu)
    }

    pub fn polarization_quark_l_thermal_part_zero_temp_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, p, mq, mu)
    }

    pub fn polarization_quark_t_thermal_part_zero_temp_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_zero_temp_landau_i(q, om, p, mq, mu)
    }

    pub fn polarization_glue_l_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> C {
        inlines::polarization_glue_l_thermal_part_landau_w_method(om, p, m, beta, integral)
    }

    pub fn polarization_glue_t_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> C {
        inlines::polarization_glue_t_thermal_part_landau_w_method(om, p, m, beta, integral)
    }

    pub fn polarization_quark_l_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
        integral: Integral,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_landau_w_method(om, p, mq, beta, mu, integral)
    }

    pub fn polarization_quark_t_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
        integral: Integral,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_landau_w_method(om, p, mq, beta, mu, integral)
    }

    pub fn polarization_quark_l_thermal_part_zero_temp_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
        integral: Integral,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_zero_temp_landau_w_method(
            om, p, mq, mu, integral,
        )
    }

    pub fn polarization_quark_t_thermal_part_zero_temp_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
        integral: Integral,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_zero_temp_landau_w_method(
            om, p, mq, mu, integral,
        )
    }

    pub fn polarization_glue_l_thermal_part_landau<T: Num>(om: T, p: R, m: R, beta: R) -> C {
        inlines::polarization_glue_l_thermal_part_landau(om, p, m, beta)
    }

    pub fn polarization_glue_t_thermal_part_landau<T: Num>(om: T, p: R, m: R, beta: R) -> C {
        inlines::polarization_glue_t_thermal_part_landau(om, p, m, beta)
    }

    pub fn polarization_quark_l_thermal_part_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_landau(om, p, mq, beta, mu)
    }

    pub fn polarization_quark_t_thermal_part_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_landau(om, p, mq, beta, mu)
    }

    pub fn polarization_quark_l_thermal_part_zero_temp_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_zero_temp_landau(om, p, mq, mu)
    }

    pub fn polarization_quark_t_thermal_part_zero_temp_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_zero_temp_landau(om, p, mq, mu)
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

    pub use super::zero_matsubara::ffi::*;
    pub use super::zero_momentum::ffi::*;

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__polarization_glue_l_thermal_part_landau_i(
        q: R,
        om: R,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::polarization_glue_l_thermal_part_landau_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__polarization_glue_t_thermal_part_landau_i(
        q: R,
        om: R,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::polarization_glue_t_thermal_part_landau_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__polarization_quark_l_thermal_part_landau_i(
        q: R,
        om: R,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_landau_i(q, om, p, mq, beta, mu)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__polarization_quark_t_thermal_part_landau_i(
        q: R,
        om: R,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_landau_i(q, om, p, mq, beta, mu)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__polarization_quark_l_thermal_part_zero_temp_landau_i(
        q: R,
        om: R,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, p, mq, mu)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__gluon__polarization_quark_t_thermal_part_zero_temp_landau_i(
        q: R,
        om: R,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        inlines::polarization_quark_t_thermal_part_zero_temp_landau_i(q, om, p, mq, mu)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::common::inlines::energy;
    use crate::common::thermal::inlines::{
        fermi_distribution_double, fermi_distribution_double_zero_temp, fermi_momentum,
    };
    use crate::consts::{get_default_integration_method, nc};
    use crate::low_level::oneloop::thermal::*;
    use crate::{Num, C, I, R};
    use peroxide::numerical::integral::{complex_integrate, Integral};
    use std::f64::consts::PI;

    const PI2: R = PI * PI;

    #[inline(always)]
    pub fn polarization_glue_l_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        let s = om * om + p * p;
        let s2 = s * s;
        let s_inv = s.inv();
        let m2 = m * m;
        let m4 = m2 * m2;
        let a = s + m2 * 2.;
        let b = s + m2;
        let b2 = b * b;

        let j_m_l_p_i = s_inv * (om * om * j_m_t_i(q, m, beta) + p * p * j_m_l_i(q, m, beta));
        let j_0_l_p_i = s_inv * (om * om * j_0_t_i(q, beta) + p * p * j_0_l_i(q, beta));
        let d_j_m_l_p_i = s_inv * (om * om * d_j_m_t_i(q, m, beta) + p * p * d_j_m_l_i(q, m, beta));
        let d2_j_m_l_p_i =
            s_inv * (om * om * d2_j_m_t_i(q, m, beta) + p * p * d2_j_m_l_i(q, m, beta));

        let t1 = (s2 * 3. / (m4 * 2.) - 1.) * i_0_0_l_i(q, om, p, beta);
        let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m4) + 4.)
            * i_m_m_l_i_same_mass(q, om, p, m, beta);
        let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_l_i(q, om, p, m, beta);
        let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, om, p, m, beta);
        let t5 = -b * s * 2. / m2 * i_m_0_i(q, om, p, m, beta);
        let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
        let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
        let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_l_i_same_mass(q, om, p, m, beta);
        let t9 = b2 / m2 * d_i_m_0_l_i(q, om, p, m, beta);
        let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, om, p, m, beta);
        let t11 = b2 * d_i_m_0_i(q, om, p, m, beta);
        let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

        let t13 = -m4 * d2_j_m_i(q, m, beta);
        let t14 = (j_m_l_p_i - j_0_l_p_i) / m2;
        let t15 = -d_j_m_l_p_i + d2_j_m_l_p_i * m2 / 2.;

        (t14 + t15 + t13) + t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11
    }

    #[inline(always)]
    pub fn polarization_glue_t_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        let s = om * om + p * p;
        let s2 = s * s;
        let m2 = m * m;
        let m4 = m2 * m2;
        let a = s + m2 * 2.;
        let b = s + m2;
        let b2 = b * b;

        let t1 = (s2 * 3. / (m4 * 2.) - 1.) * i_0_0_t_i(q, om, p, beta);
        let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m4) + 4.)
            * i_m_m_t_i_same_mass(q, om, p, m, beta);
        let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_t_i(q, om, p, m, beta);
        let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, om, p, m, beta);
        let t5 = -b * s * 2. / m2 * i_m_0_i(q, om, p, m, beta);
        let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
        let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
        let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_t_i_same_mass(q, om, p, m, beta);
        let t9 = b2 / m2 * d_i_m_0_t_i(q, om, p, m, beta);
        let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, om, p, m, beta);
        let t11 = b2 * d_i_m_0_i(q, om, p, m, beta);
        let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

        let t13 = -m4 * d2_j_m_i(q, m, beta);
        let t14 = (j_m_t_i(q, m, beta) - j_0_t_i(q, beta)) / m2;
        let t15 = -d_j_m_t_i(q, m, beta) + d2_j_m_t_i(q, m, beta) * m2 / 2.;

        t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11 + t14 + t15 + t13
    }

    #[inline(always)]
    pub fn polarization_quark_l_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        let p2 = p * p;
        let s = om * om + p2;
        let en = energy(q, mq);
        let tl = tlog_same_mass(q, om, p, en);
        let tl_opp = tlog_same_mass(q, -om, p, en);
        let t1 = (s - en * en * 4. + om * 4. * I * en) * tl;
        let t1_opp = (s - en * en * 4. - om * 4. * I * en) * tl_opp;
        -s * fermi_distribution_double(en, beta, mu) * q * q / (p2 * en * (nc() as R) * 2. * PI2)
            * (1. - (t1 + t1_opp) / (8. * q * p))
    }

    #[inline(always)]
    pub fn polarization_quark_t_thermal_part_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        let p2 = p * p;
        let om2 = om * om;
        let s = om2 + p2;
        let q2 = q * q;
        let en = energy(q, mq);
        let a = -om2 + p2;
        let b = a * s + 4. * q2 * p2 + om2 * 4. * en * en;
        let c = om * s * 4. * I * en;
        let tl = tlog_same_mass(q, om, p, en);
        let tl_opp = tlog_same_mass(q, -om, p, en);
        let t1 = (b - c) * tl;
        let t1_opp = (b + c) * tl_opp;
        -fermi_distribution_double(en, beta, mu) * q2 / (p2 * en * (nc() as R) * 4. * PI2)
            * (a - (t1 + t1_opp) / (8. * q * p))
    }

    #[inline(always)]
    pub fn polarization_quark_l_thermal_part_zero_temp_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        let p2 = p * p;
        let s = om * om + p2;
        let en = energy(q, mq);
        let tl = tlog_same_mass(q, om, p, en);
        let tl_opp = tlog_same_mass(q, -om, p, en);
        let t1 = (s - en * en * 4. + om * 4. * I * en) * tl;
        let t1_opp = (s - en * en * 4. - om * 4. * I * en) * tl_opp;
        -s * fermi_distribution_double_zero_temp(en, mu) * q * q
            / (p2 * en * (nc() as R) * 2. * PI2)
            * (1. - (t1 + t1_opp) / (8. * q * p))
    }

    #[inline(always)]
    pub fn polarization_quark_t_thermal_part_zero_temp_landau_i<T: Num>(
        q: R,
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        let p2 = p * p;
        let om2 = om * om;
        let s = om2 + p2;
        let q2 = q * q;
        let en = energy(q, mq);
        let a = -om2 + p2;
        let b = a * s + 4. * q2 * p2 + om2 * 4. * en * en;
        let c = om * s * 4. * I * en;
        let tl = tlog_same_mass(q, om, p, en);
        let tl_opp = tlog_same_mass(q, -om, p, en);
        let t1 = (b - c) * tl;
        let t1_opp = (b + c) * tl_opp;
        -fermi_distribution_double_zero_temp(en, mu) * q2 / (p2 * en * (nc() as R) * 4. * PI2)
            * (a - (t1 + t1_opp) / (8. * q * p))
    }

    #[inline(always)]
    pub fn polarization_glue_l_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> C {
        complex_integrate(
            |t| polarization_glue_l_thermal_part_landau_i((1. - t) / t, om, p, m, beta) / (t * t),
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn polarization_glue_t_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> C {
        complex_integrate(
            |t| polarization_glue_t_thermal_part_landau_i((1. - t) / t, om, p, m, beta) / (t * t),
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn polarization_quark_l_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
        integral: Integral,
    ) -> C {
        complex_integrate(
            |t| {
                polarization_quark_l_thermal_part_landau_i((1. - t) / t, om, p, mq, beta, mu)
                    / (t * t)
            },
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn polarization_quark_t_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
        integral: Integral,
    ) -> C {
        complex_integrate(
            |t| {
                polarization_quark_t_thermal_part_landau_i((1. - t) / t, om, p, mq, beta, mu)
                    / (t * t)
            },
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn polarization_quark_l_thermal_part_zero_temp_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
        integral: Integral,
    ) -> C {
        if mu.abs() < mq.abs() {
            return C::new(0., 0.);
        };
        complex_integrate(
            |q| polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, p, mq, mu),
            (0., fermi_momentum(mq, mu)),
            integral,
        )
    }

    #[inline(always)]
    pub fn polarization_quark_t_thermal_part_zero_temp_landau_w_method<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
        integral: Integral,
    ) -> C {
        if mu.abs() < mq.abs() {
            return C::new(0., 0.);
        };
        complex_integrate(
            |q| polarization_quark_t_thermal_part_zero_temp_landau_i(q, om, p, mq, mu),
            (0., fermi_momentum(mq, mu)),
            integral,
        )
    }

    #[inline(always)]
    pub fn polarization_glue_l_thermal_part_landau<T: Num>(om: T, p: R, m: R, beta: R) -> C {
        polarization_glue_l_thermal_part_landau_w_method(
            om,
            p,
            m,
            beta,
            get_default_integration_method(),
        )
    }

    #[inline(always)]
    pub fn polarization_glue_t_thermal_part_landau<T: Num>(om: T, p: R, m: R, beta: R) -> C {
        polarization_glue_t_thermal_part_landau_w_method(
            om,
            p,
            m,
            beta,
            get_default_integration_method(),
        )
    }

    #[inline(always)]
    pub fn polarization_quark_l_thermal_part_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        polarization_quark_l_thermal_part_landau_w_method(
            om,
            p,
            mq,
            beta,
            mu,
            get_default_integration_method(),
        )
    }

    #[inline(always)]
    pub fn polarization_quark_t_thermal_part_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        beta: R,
        mu: R,
    ) -> C {
        polarization_quark_t_thermal_part_landau_w_method(
            om,
            p,
            mq,
            beta,
            mu,
            get_default_integration_method(),
        )
    }

    #[inline(always)]
    pub fn polarization_quark_l_thermal_part_zero_temp_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        polarization_quark_l_thermal_part_zero_temp_landau_w_method(
            om,
            p,
            mq,
            mu,
            get_default_integration_method(),
        )
    }

    #[inline(always)]
    pub fn polarization_quark_t_thermal_part_zero_temp_landau<T: Num>(
        om: T,
        p: R,
        mq: R,
        mu: R,
    ) -> C {
        polarization_quark_t_thermal_part_zero_temp_landau_w_method(
            om,
            p,
            mq,
            mu,
            get_default_integration_method(),
        )
    }
}

pub mod zero_matsubara {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::Integral;
        use crate::{C, R};

        pub fn polarization_glue_l_thermal_part_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::polarization_glue_l_thermal_part_landau_i(q, p, m, beta)
        }

        pub fn polarization_glue_t_thermal_part_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::polarization_glue_t_thermal_part_landau_i(q, p, m, beta)
        }

        pub fn polarization_quark_l_thermal_part_landau_i(q: R, p: R, mq: R, beta: R, mu: R) -> R {
            inlines::polarization_quark_l_thermal_part_landau_i(q, p, mq, beta, mu)
        }

        pub fn polarization_quark_t_thermal_part_landau_i(q: R, p: R, mq: R, beta: R, mu: R) -> R {
            inlines::polarization_quark_t_thermal_part_landau_i(q, p, mq, beta, mu)
        }

        pub fn polarization_quark_l_thermal_part_zero_temp_landau_i(q: R, p: R, mq: R, mu: R) -> R {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau_i(q, p, mq, mu)
        }

        pub fn polarization_quark_t_thermal_part_zero_temp_landau_i(q: R, p: R, mq: R, mu: R) -> R {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau_i(q, p, mq, mu)
        }

        pub fn polarization_glue_l_thermal_part_landau_w_method(
            p: R,
            m: R,
            beta: R,
            integral: Integral,
        ) -> R {
            inlines::polarization_glue_l_thermal_part_landau_w_method(p, m, beta, integral)
        }

        pub fn polarization_glue_t_thermal_part_landau_w_method(
            p: R,
            m: R,
            beta: R,
            integral: Integral,
        ) -> R {
            inlines::polarization_glue_t_thermal_part_landau_w_method(p, m, beta, integral)
        }

        pub fn polarization_quark_l_thermal_part_landau_w_method(
            p: R,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> R {
            inlines::polarization_quark_l_thermal_part_landau_w_method(p, mq, beta, mu, integral)
        }

        pub fn polarization_quark_t_thermal_part_landau_w_method(
            p: R,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> R {
            inlines::polarization_quark_t_thermal_part_landau_w_method(p, mq, beta, mu, integral)
        }

        pub fn polarization_quark_l_thermal_part_zero_temp_landau_w_method(
            p: R,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> R {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau_w_method(
                p, mq, mu, integral,
            )
        }

        pub fn polarization_quark_t_thermal_part_zero_temp_landau_w_method(
            p: R,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> R {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau_w_method(
                p, mq, mu, integral,
            )
        }

        pub fn polarization_glue_l_thermal_part_landau(p: R, m: R, beta: R) -> R {
            inlines::polarization_glue_l_thermal_part_landau(p, m, beta)
        }

        pub fn polarization_glue_t_thermal_part_landau(p: R, m: R, beta: R) -> R {
            inlines::polarization_glue_t_thermal_part_landau(p, m, beta)
        }

        pub fn polarization_quark_l_thermal_part_landau(p: R, mq: R, beta: R, mu: R) -> R {
            inlines::polarization_quark_l_thermal_part_landau(p, mq, beta, mu)
        }

        pub fn polarization_quark_t_thermal_part_landau(p: R, mq: R, beta: R, mu: R) -> R {
            inlines::polarization_quark_t_thermal_part_landau(p, mq, beta, mu)
        }

        pub fn polarization_quark_l_thermal_part_zero_temp_landau(p: R, mq: R, mu: R) -> R {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau(p, mq, mu)
        }

        pub fn polarization_quark_t_thermal_part_zero_temp_landau(p: R, mq: R, mu: R) -> R {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau(p, mq, mu)
        }
    }

    pub(crate) mod ffi {
        use super::inlines;
        use crate::{C, R};

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__gluon__polarization_glue_l_thermal_part_landau_i(
            q: R,
            om: R,
            m: R,
            beta: R,
        ) -> C {
            inlines::polarization_glue_l_thermal_part_landau_i(q, om, m, beta)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__gluon__polarization_glue_t_thermal_part_landau_i(
            q: R,
            om: R,
            m: R,
            beta: R,
        ) -> C {
            inlines::polarization_glue_t_thermal_part_landau_i(q, om, m, beta)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__gluon__polarization_quark_l_thermal_part_landau_i(
            q: R,
            p: R,
            mq: R,
            beta: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_l_thermal_part_landau_i(q, p, mq, beta, mu)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__gluon__polarization_quark_t_thermal_part_landau_i(
            q: R,
            p: R,
            mq: R,
            beta: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_t_thermal_part_landau_i(q, p, mq, beta, mu)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__gluon__polarization_quark_l_thermal_part_zero_temp_landau_i(
            q: R,
            p: R,
            mq: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau_i(q, p, mq, mu)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__gluon__polarization_quark_t_thermal_part_zero_temp_landau_i(
            q: R,
            p: R,
            mq: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau_i(q, p, mq, mu)
        }
    }

    mod inlines {
        use crate::common::inlines::energy;
        use crate::common::thermal::inlines::{
            fermi_distribution_double, fermi_distribution_double_zero_temp, fermi_momentum,
        };
        use crate::consts::{
            get_default_integration_method, matsubara_reg, max_iter_integral, nc, tol_integral,
        };
        use crate::low_level::oneloop::thermal::zero_matsubara::*;
        use crate::low_level::oneloop::thermal::{
            d2_j_m_i, d2_j_m_l_i, d2_j_m_t_i, d_j_m_i, d_j_m_l_i, d_j_m_t_i, j_0_i, j_0_l_i,
            j_0_t_i, j_m_i, j_m_l_i, j_m_t_i,
        };
        use crate::{C, R};
        use peroxide::numerical::integral::{integrate, Integral};
        use std::f64::consts::PI;

        const DEFAULT_TOL_INT_REPL: R = 1E-8;
        const PI2: R = PI * PI;

        #[inline(always)]
        pub fn polarization_glue_l_thermal_part_landau_i(q: R, p: R, m: R, beta: R) -> C {
            let s = p * p;
            let s2 = s * s;
            let m2 = m * m;
            let m4 = m2 * m2;
            let a = m2.mul_add(2., s);
            let b = s + m2;

            let j_m_l_p_i = j_m_l_i(q, m, beta);
            let j_0_l_p_i = j_0_l_i(q, beta);
            let d_j_m_l_p_i = d_j_m_l_i(q, m, beta);
            let d2_j_m_l_p_i = d2_j_m_l_i(q, m, beta);

            let t1 = (s2 * 3. / (m4 * 2.) - 1.) * i_0_0_l_i(q, p, beta);
            let t2 = (4.0f64.mul_add(m4, s2.mul_add(3., s * m2 * 8.)) / (2. * m4) + 4.)
                * i_m_m_l_i_same_mass(q, p, m, beta);
            let t3 = -(s2.mul_add(3., s * m2 * 4.) + m4) / m4 * i_m_0_l_i(q, p, m, beta);
            let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, p, m, beta);
            let t5 = -b * s * 2. / m2 * i_m_0_i(q, p, m, beta);
            let t6 = -s.mul_add(2., m2 * 3.) / m2 * j_m_i(q, m, beta);
            let t7 = s.mul_add(2., m2) / m2 * j_0_i(q, beta);
            let t8 = -m2.mul_add(8., a * a / m2) * d_i_m_m_l_i_same_mass(q, p, m, beta);
            let t9 = b * b / m2 * d_i_m_0_l_i(q, p, m, beta);
            let t10 = -s * 2. * m2.mul_add(4., s) * d_i_m_m_i_same_mass(q, p, m, beta);
            let t11 = b * b * d_i_m_0_i(q, p, m, beta);
            let t12 = m2.mul_add(3., s) * d_j_m_i(q, m, beta);

            let t13 = -m4 * d2_j_m_i(q, m, beta);
            let t14 = (j_m_l_p_i - j_0_l_p_i) / m2;
            let t15 = -d_j_m_l_p_i + d2_j_m_l_p_i * m2 / 2.;

            (t14 + t15 + t13) + t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11
        }

        #[inline(always)]
        pub fn polarization_glue_t_thermal_part_landau_i(q: R, p: R, m: R, beta: R) -> C {
            let s = p * p;
            let s2 = s * s;
            let m2 = m * m;
            let m4 = m2 * m2;
            let a = m2.mul_add(2., s);
            let b = s + m2;

            let t1 = (s2 * 3. / (m4 * 2.) - 1.) * i_0_0_t_i(q, p, beta);
            let t2 = (4.0f64.mul_add(m4, s2.mul_add(3., s * m2 * 8.)) / (2. * m4) + 4.)
                * i_m_m_t_i_same_mass(q, p, m, beta);
            let t3 = -(s2.mul_add(3., s * m2 * 4.) + m4) / m4 * i_m_0_t_i(q, p, m, beta);
            let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, p, m, beta);
            let t5 = -b * s * 2. / m2 * i_m_0_i(q, p, m, beta);
            let t6 = -s.mul_add(2., m2 * 3.) / m2 * j_m_i(q, m, beta);
            let t7 = s.mul_add(2., m2) / m2 * j_0_i(q, beta);
            let t8 = -m2.mul_add(8., a * a / m2) * d_i_m_m_t_i_same_mass(q, p, m, beta);
            let t9 = b * b / m2 * d_i_m_0_t_i(q, p, m, beta);
            let t10 = -s * 2. * m2.mul_add(4., s) * d_i_m_m_i_same_mass(q, p, m, beta);
            let t11 = b * b * d_i_m_0_i(q, p, m, beta);
            let t12 = m2.mul_add(3., s) * d_j_m_i(q, m, beta);

            let t13 = -m4 * d2_j_m_i(q, m, beta);
            let t14 = (j_m_t_i(q, m, beta) - j_0_t_i(q, beta)) / m2;
            let t15 = -d_j_m_t_i(q, m, beta) + d2_j_m_t_i(q, m, beta) * m2 / 2.;

            t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11 + t14 + t15 + t13
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_landau_i(q: R, p: R, mq: R, beta: R, mu: R) -> R {
            let s = p * p;
            let en = energy(q, mq);
            let t1 = (en * en).mul_add(-4., s) / (4. * q * p) * tlog_same_mass(q, p);
            -fermi_distribution_double(en, beta, mu) * q * q / (en * (nc() as R) * 2. * PI2)
                * (1. - t1)
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_landau_i(q: R, p: R, mq: R, beta: R, mu: R) -> R {
            let p2 = p * p;
            let q2 = q * q;
            let en = energy(q, mq);
            let t1 = 4.0f64.mul_add(q2, p2) * tlog_same_mass(q, p);
            -fermi_distribution_double(en, beta, mu) * q2 / (en * (nc() as R) * 4. * PI2)
                * (1. - t1 / (4. * q * p))
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_zero_temp_landau_i(q: R, p: R, mq: R, mu: R) -> R {
            let s = p * p;
            let en = energy(q, mq);
            let t1 = (en * en).mul_add(-4., s) / (4. * q * p) * tlog_same_mass(q, p);
            -fermi_distribution_double_zero_temp(en, mu) * q * q / (en * (nc() as R) * 2. * PI2)
                * (1. - t1)
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_zero_temp_landau_i(q: R, p: R, mq: R, mu: R) -> R {
            let p2 = p * p;
            let q2 = q * q;
            let en = energy(q, mq);
            let t1 = 4.0f64.mul_add(q2, p2) * tlog_same_mass(q, p);
            -fermi_distribution_double_zero_temp(en, mu) * q2 / (en * (nc() as R) * 4. * PI2)
                * (1. - t1 / (4. * q * p))
        }

        #[inline(always)]
        pub fn polarization_glue_l_thermal_part_landau_w_method(
            p: R,
            m: R,
            beta: R,
            integral: Integral,
        ) -> R {
            // You may be thinking: did you spend countless ours of your life
            // computing and coding the Matsubara limits to then end up using
            // the ordinary thermal functions with omega=epsilon?! The answer
            // is... yeah, 'cos apparently there's a perfectly-well-behaved
            // logarithmic singularity in the integrands in the Matsubara limit
            // which however fucks around with numerical integration and yields
            // plots which look very much like the ECG of my grandma with
            // arrythmia.
            integrate(
                |t| {
                    super::super::polarization_glue_l_thermal_part_landau_i(
                        (1. - t) / t,
                        matsubara_reg(),
                        p,
                        m,
                        beta,
                    )
                    .re / (t * t)
                },
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_glue_t_thermal_part_landau_w_method(
            p: R,
            m: R,
            beta: R,
            integral: Integral,
        ) -> R {
            integrate(
                |t| {
                    super::super::polarization_glue_t_thermal_part_landau_i(
                        (1. - t) / t,
                        matsubara_reg(),
                        p,
                        m,
                        beta,
                    )
                    .re / (t * t)
                },
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_landau_w_method(
            p: R,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> R {
            integrate(
                |t| {
                    super::super::polarization_quark_l_thermal_part_landau_i(
                        (1. - t) / t,
                        matsubara_reg(),
                        p,
                        mq,
                        beta,
                        mu,
                    )
                    .re / (t * t)
                },
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_landau_w_method(
            p: R,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> R {
            integrate(
                |t| {
                    super::super::polarization_quark_t_thermal_part_landau_i(
                        (1. - t) / t,
                        matsubara_reg(),
                        p,
                        mq,
                        beta,
                        mu,
                    )
                    .re / (t * t)
                },
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_zero_temp_landau_w_method(
            p: R,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> R {
            if mu.abs() < mq.abs() {
                return 0.;
            };
            integrate(
                |q| {
                    super::super::polarization_quark_l_thermal_part_zero_temp_landau_i(
                        q,
                        matsubara_reg(),
                        p,
                        mq,
                        mu,
                    )
                    .re
                },
                (0., fermi_momentum(mq, mu)),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_zero_temp_landau_w_method(
            p: R,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> R {
            if mu.abs() < mq.abs() {
                return 0.;
            };
            integrate(
                |q| {
                    super::super::polarization_quark_t_thermal_part_zero_temp_landau_i(
                        q,
                        matsubara_reg(),
                        p,
                        mq,
                        mu,
                    )
                    .re
                },
                (0., fermi_momentum(mq, mu)),
                integral,
            )
        }

        #[inline(always)]
        // To make things worse, if EPS (the regulator for the matsubara frequency)
        // is too small, the integral's convergence becomes very slow. So, unless the
        // user customized the integration method, we increase the tolerance of the
        // default integration method to DEFAULT_TOL_INT_REPL. This does not increase
        // the error meaningfully, since the error due to EPS being non zero is larger
        // by at least two orders of magnitude (given the defaults at the time of writing).
        pub fn polarization_glue_l_thermal_part_landau(p: R, m: R, beta: R) -> R {
            let (mut integral, dti) = (get_default_integration_method(), tol_integral());
            if let crate::Integral::G7K15(ti, mi) = integral {
                #[allow(clippy::float_cmp)]
                if ti == dti && mi == max_iter_integral() {
                    integral = crate::Integral::G7K15(DEFAULT_TOL_INT_REPL, mi);
                }
            }
            polarization_glue_l_thermal_part_landau_w_method(p, m, beta, integral)
        }

        #[inline(always)]
        pub fn polarization_glue_t_thermal_part_landau(p: R, m: R, beta: R) -> R {
            let (mut integral, dti) = (get_default_integration_method(), tol_integral());
            if let crate::Integral::G7K15(ti, mi) = integral {
                #[allow(clippy::float_cmp)]
                if ti == dti && mi == max_iter_integral() {
                    integral = crate::Integral::G7K15(DEFAULT_TOL_INT_REPL, mi);
                }
            }
            polarization_glue_t_thermal_part_landau_w_method(p, m, beta, integral)
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_landau(p: R, mq: R, beta: R, mu: R) -> R {
            let (mut integral, dti) = (get_default_integration_method(), tol_integral());
            if let crate::Integral::G7K15(ti, mi) = integral {
                #[allow(clippy::float_cmp)]
                if ti == dti && mi == max_iter_integral() {
                    integral = crate::Integral::G7K15(DEFAULT_TOL_INT_REPL, mi);
                }
            }
            polarization_quark_l_thermal_part_landau_w_method(p, mq, beta, mu, integral)
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_landau(p: R, mq: R, beta: R, mu: R) -> R {
            let (mut integral, dti) = (get_default_integration_method(), tol_integral());
            if let crate::Integral::G7K15(ti, mi) = integral {
                #[allow(clippy::float_cmp)]
                if ti == dti && mi == max_iter_integral() {
                    integral = crate::Integral::G7K15(DEFAULT_TOL_INT_REPL, mi);
                }
            }
            polarization_quark_t_thermal_part_landau_w_method(p, mq, beta, mu, integral)
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_zero_temp_landau(p: R, mq: R, mu: R) -> R {
            let (mut integral, dti) = (get_default_integration_method(), tol_integral());
            if let crate::Integral::G7K15(ti, mi) = integral {
                #[allow(clippy::float_cmp)]
                if ti == dti && mi == max_iter_integral() {
                    integral = crate::Integral::G7K15(DEFAULT_TOL_INT_REPL, mi);
                }
            }
            polarization_quark_l_thermal_part_zero_temp_landau_w_method(p, mq, mu, integral)
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_zero_temp_landau(p: R, mq: R, mu: R) -> R {
            let (mut integral, dti) = (get_default_integration_method(), tol_integral());
            if let crate::Integral::G7K15(ti, mi) = integral {
                #[allow(clippy::float_cmp)]
                if ti == dti && mi == max_iter_integral() {
                    integral = crate::Integral::G7K15(DEFAULT_TOL_INT_REPL, mi);
                }
            }
            polarization_quark_t_thermal_part_zero_temp_landau_w_method(p, mq, mu, integral)
        }
    }
}

pub mod zero_momentum {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::Integral;
        use crate::{Num, C, R};

        pub fn polarization_glue_l_thermal_part_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            inlines::polarization_glue_l_thermal_part_landau_i(q, om, m, beta)
        }

        pub fn polarization_glue_t_thermal_part_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            inlines::polarization_glue_t_thermal_part_landau_i(q, om, m, beta)
        }

        pub fn polarization_quark_l_thermal_part_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            beta: R,
            mu: R,
        ) -> T {
            inlines::polarization_quark_l_thermal_part_landau_i(q, om, mq, beta, mu)
        }

        pub fn polarization_quark_t_thermal_part_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            beta: R,
            mu: R,
        ) -> T {
            inlines::polarization_quark_t_thermal_part_landau_i(q, om, mq, beta, mu)
        }

        pub fn polarization_quark_l_thermal_part_zero_temp_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            mu: R,
        ) -> T {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, mq, mu)
        }

        pub fn polarization_quark_t_thermal_part_zero_temp_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            mu: R,
        ) -> T {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau_i(q, om, mq, mu)
        }

        pub fn polarization_glue_l_thermal_part_landau_w_method<T: Num>(
            om: T,
            m: R,
            beta: R,
            integral: Integral,
        ) -> C {
            inlines::polarization_glue_l_thermal_part_landau_w_method(om, m, beta, integral)
        }

        pub fn polarization_glue_t_thermal_part_landau_w_method<T: Num>(
            om: T,
            m: R,
            beta: R,
            integral: Integral,
        ) -> C {
            inlines::polarization_glue_t_thermal_part_landau_w_method(om, m, beta, integral)
        }

        pub fn polarization_quark_l_thermal_part_landau_w_method<T: Num>(
            om: T,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> C {
            inlines::polarization_quark_l_thermal_part_landau_w_method(om, mq, beta, mu, integral)
        }

        pub fn polarization_quark_t_thermal_part_landau_w_method<T: Num>(
            om: T,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> C {
            inlines::polarization_quark_t_thermal_part_landau_w_method(om, mq, beta, mu, integral)
        }

        pub fn polarization_quark_l_thermal_part_zero_temp_landau_w_method<T: Num>(
            om: T,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> C {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau_w_method(
                om, mq, mu, integral,
            )
        }

        pub fn polarization_quark_t_thermal_part_zero_temp_landau_w_method<T: Num>(
            om: T,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> C {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau_w_method(
                om, mq, mu, integral,
            )
        }

        pub fn polarization_glue_l_thermal_part_landau<T: Num>(om: T, m: R, beta: R) -> C {
            inlines::polarization_glue_l_thermal_part_landau(om, m, beta)
        }

        pub fn polarization_glue_t_thermal_part_landau<T: Num>(om: T, m: R, beta: R) -> C {
            inlines::polarization_glue_t_thermal_part_landau(om, m, beta)
        }

        pub fn polarization_quark_l_thermal_part_landau<T: Num>(om: T, mq: R, beta: R, mu: R) -> C {
            inlines::polarization_quark_l_thermal_part_landau(om, mq, beta, mu)
        }

        pub fn polarization_quark_t_thermal_part_landau<T: Num>(om: T, mq: R, beta: R, mu: R) -> C {
            inlines::polarization_quark_t_thermal_part_landau(om, mq, beta, mu)
        }

        pub fn polarization_quark_l_thermal_part_zero_temp_landau<T: Num>(
            om: T,
            mq: R,
            mu: R,
        ) -> C {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau(om, mq, mu)
        }

        pub fn polarization_quark_t_thermal_part_zero_temp_landau<T: Num>(
            om: T,
            mq: R,
            mu: R,
        ) -> C {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau(om, mq, mu)
        }
    }

    pub(crate) mod ffi {
        use super::inlines;
        use crate::{C, R};

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__gluon__polarization_glue_l_thermal_part_landau_i(
            q: R,
            om: R,
            m: R,
            beta: R,
        ) -> C {
            inlines::polarization_glue_l_thermal_part_landau_i(q, om, m, beta)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__gluon__polarization_glue_t_thermal_part_landau_i(
            q: R,
            om: R,
            m: R,
            beta: R,
        ) -> C {
            inlines::polarization_glue_t_thermal_part_landau_i(q, om, m, beta)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__gluon__polarization_quark_l_thermal_part_landau_i(
            q: R,
            om: R,
            mq: R,
            beta: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_l_thermal_part_landau_i(q, om, mq, beta, mu)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__gluon__polarization_quark_t_thermal_part_landau_i(
            q: R,
            om: R,
            mq: R,
            beta: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_t_thermal_part_landau_i(q, om, mq, beta, mu)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__gluon__polarization_quark_l_thermal_part_zero_temp_landau_i(
            q: R,
            om: R,
            mq: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, mq, mu)
        }

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__gluon__polarization_quark_t_thermal_part_zero_temp_landau_i(
            q: R,
            om: R,
            mq: R,
            mu: R,
        ) -> R {
            inlines::polarization_quark_t_thermal_part_zero_temp_landau_i(q, om, mq, mu)
        }
    }

    mod inlines {
        use crate::common::inlines::energy;
        use crate::common::thermal::inlines::{
            fermi_distribution_double, fermi_distribution_double_zero_temp,
        };
        use crate::consts::{get_default_integration_method, nc};
        use crate::low_level::oneloop::thermal::zero_momentum::*;
        use crate::low_level::oneloop::thermal::{
            d2_j_m_i, d2_j_m_t_i, d_j_m_i, d_j_m_t_i, j_0_i, j_0_t_i, j_m_i, j_m_t_i,
        };
        use crate::{Num, C, R};
        use peroxide::numerical::integral::{complex_integrate, Integral};
        use std::f64::consts::PI;

        const PI2: R = PI * PI;

        #[inline(always)]
        pub fn polarization_glue_l_thermal_part_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            let s = om * om;
            let s2 = s * s;
            let m2 = m * m;
            let m4 = m2 * m2;
            let a = s + m2 * 2.;
            let b = s + m2;

            let j_m_l_p_i = j_m_t_i(q, m, beta);
            let j_0_l_p_i = j_0_t_i(q, beta);
            let d_j_m_l_p_i = d_j_m_t_i(q, m, beta);
            let d2_j_m_l_p_i = d2_j_m_t_i(q, m, beta);

            let t1 = (s2 * 3. / (m4 * 2.) - 1.) * i_0_0_l_i(q, om, beta);
            let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m4) + 4.)
                * i_m_m_l_i_same_mass(q, om, m, beta);
            let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_l_i(q, om, m, beta);
            let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, om, m, beta);
            let t5 = -b * s * 2. / m2 * i_m_0_i(q, om, m, beta);
            let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
            let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
            let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_l_i_same_mass(q, om, m, beta);
            let t9 = b * b / m2 * d_i_m_0_l_i(q, om, m, beta);
            let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, om, m, beta);
            let t11 = b * b * d_i_m_0_i(q, om, m, beta);
            let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

            let t13 = -m4 * d2_j_m_i(q, m, beta);
            let t14 = (j_m_l_p_i - j_0_l_p_i) / m2;
            let t15 = -d_j_m_l_p_i + d2_j_m_l_p_i * m2 / 2.;

            t12 + (t14 + t15 + t13) + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11
        }

        #[inline(always)]
        pub fn polarization_glue_t_thermal_part_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            polarization_glue_l_thermal_part_landau_i(q, om, m, beta)
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            beta: R,
            mu: R,
        ) -> T {
            let en = energy(q, mq);
            let en2 = en * en;
            let q2 = q * q;

            ((om * om + 4. * en2) * en * 3. * PI2 * (nc() as R)).inv()
                * 2.
                * fermi_distribution_double(en, beta, mu)
                * q2
                * 3.0f64.mul_add(-en2, q2)
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            beta: R,
            mu: R,
        ) -> T {
            polarization_quark_l_thermal_part_landau_i(q, om, mq, beta, mu)
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_zero_temp_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            mu: R,
        ) -> T {
            let en = energy(q, mq);
            let en2 = en * en;
            let q2 = q * q;

            ((om * om + 4. * en2) * en * 3. * PI2 * (nc() as R)).inv()
                * 2.
                * fermi_distribution_double_zero_temp(en, mu)
                * q2
                * 3.0f64.mul_add(-en2, q2)
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_zero_temp_landau_i<T: Num>(
            q: R,
            om: T,
            mq: R,
            mu: R,
        ) -> T {
            polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, mq, mu)
        }

        #[inline(always)]
        pub fn polarization_glue_l_thermal_part_landau_w_method<T: Num>(
            om: T,
            m: R,
            beta: R,
            integral: Integral,
        ) -> C {
            complex_integrate(
                |t| polarization_glue_l_thermal_part_landau_i((1. - t) / t, om, m, beta) / (t * t),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_glue_t_thermal_part_landau_w_method<T: Num>(
            om: T,
            m: R,
            beta: R,
            integral: Integral,
        ) -> C {
            complex_integrate(
                |t| polarization_glue_t_thermal_part_landau_i((1. - t) / t, om, m, beta) / (t * t),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_landau_w_method<T: Num>(
            om: T,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> C {
            complex_integrate(
                |t| {
                    polarization_quark_l_thermal_part_landau_i((1. - t) / t, om, mq, beta, mu)
                        .into()
                        / (t * t)
                },
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_landau_w_method<T: Num>(
            om: T,
            mq: R,
            beta: R,
            mu: R,
            integral: Integral,
        ) -> C {
            complex_integrate(
                |t| {
                    polarization_quark_t_thermal_part_landau_i((1. - t) / t, om, mq, beta, mu)
                        .into()
                        / (t * t)
                },
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_zero_temp_landau_w_method<T: Num>(
            om: T,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> C {
            if mu.abs() < mq.abs() {
                return C::new(0., 0.);
            };
            complex_integrate(
                |q| polarization_quark_l_thermal_part_zero_temp_landau_i(q, om, mq, mu).into(),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_zero_temp_landau_w_method<T: Num>(
            om: T,
            mq: R,
            mu: R,
            integral: Integral,
        ) -> C {
            if mu.abs() < mq.abs() {
                return C::new(0., 0.);
            };
            complex_integrate(
                |q| polarization_quark_t_thermal_part_zero_temp_landau_i(q, om, mq, mu).into(),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn polarization_glue_l_thermal_part_landau<T: Num>(om: T, m: R, beta: R) -> C {
            polarization_glue_l_thermal_part_landau_w_method(
                om,
                m,
                beta,
                get_default_integration_method(),
            )
        }

        #[inline(always)]
        pub fn polarization_glue_t_thermal_part_landau<T: Num>(om: T, m: R, beta: R) -> C {
            polarization_glue_t_thermal_part_landau_w_method(
                om,
                m,
                beta,
                get_default_integration_method(),
            )
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_landau<T: Num>(om: T, mq: R, beta: R, mu: R) -> C {
            polarization_quark_l_thermal_part_landau_w_method(
                om,
                mq,
                beta,
                mu,
                get_default_integration_method(),
            )
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_landau<T: Num>(om: T, mq: R, beta: R, mu: R) -> C {
            polarization_quark_t_thermal_part_landau_w_method(
                om,
                mq,
                beta,
                mu,
                get_default_integration_method(),
            )
        }

        #[inline(always)]
        pub fn polarization_quark_l_thermal_part_zero_temp_landau<T: Num>(
            om: T,
            mq: R,
            mu: R,
        ) -> C {
            polarization_quark_l_thermal_part_zero_temp_landau_w_method(
                om,
                mq,
                mu,
                get_default_integration_method(),
            )
        }

        #[inline(always)]
        pub fn polarization_quark_t_thermal_part_zero_temp_landau<T: Num>(
            om: T,
            mq: R,
            mu: R,
        ) -> C {
            polarization_quark_t_thermal_part_zero_temp_landau_w_method(
                om,
                mq,
                mu,
                get_default_integration_method(),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Num, C, I, R};

    const TOLERANCE: R = 1e-11;

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
    fn test_polarization_glue_l_thermal_part_landau_i() {
        use super::{
            ffi::oneloop__gluon__polarization_glue_l_thermal_part_landau_i,
            polarization_glue_l_thermal_part_landau_i,
        };

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.0017852067865610625 + 0. * I,
            -0.0030398086246881196 + 0. * I,
            -0.0026809417392170062 + 0. * I,
            -0.0023335168072622214 + 0. * I,
            0.0006964565160878051 + 0. * I,
            -0.8974644543820465 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    polarization_glue_l_thermal_part_landau_i(*q, *om, *p, *m, *beta),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    oneloop__gluon__polarization_glue_l_thermal_part_landau_i(
                        *q, *om, *p, *m, *beta,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_polarization_glue_t_thermal_part_landau_i() {
        use super::{
            ffi::oneloop__gluon__polarization_glue_t_thermal_part_landau_i,
            polarization_glue_t_thermal_part_landau_i,
        };

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.004364759744835432 + 0. * I,
            -0.003176656270982809 + 0. * I,
            -0.0022201869715893826 + 0. * I,
            -0.0017757726773228275 + 0. * I,
            -0.0004258083320350534 + 0. * I,
            -0.19221349678773383 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    polarization_glue_t_thermal_part_landau_i(*q, *om, *p, *m, *beta),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    oneloop__gluon__polarization_glue_t_thermal_part_landau_i(
                        *q, *om, *p, *m, *beta,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_polarization_quark_l_thermal_part_landau_i() {
        use super::{
            ffi::oneloop__gluon__polarization_quark_l_thermal_part_landau_i,
            polarization_quark_l_thermal_part_landau_i,
        };

        assert_eq!(crate::consts::nc(), 3);

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 3.2, 0.8),
            (0.35, 0.21, 2.16, 1.2, 3.2, 0.8),
            (0.35, 0.75, 2.16, 1.2, 3.2, 0.8),
            (0.35, 0.75, 1.15, 1.2, 3.2, 0.8),
            (0.35, 0.75, 1.15, 3.76, 3.2, 0.8),
            (0.35, 0.75, 1.15, 3.76, 0.19, 0.8),
            (0.35, 0.75, 1.15, 3.76, 0.19, 5.9),
        ];

        let res: [C; 7] = [
            -0.001120552759944293 + 0. * I,
            -0.0004212440449496511 + 0. * I,
            -0.0003336641405716638 + 0. * I,
            -0.0005021727489952145 + 0. * I,
            -1.1925807073568502e-07 + 0. * I,
            -0.0010669564506489019 + 0. * I,
            -0.001195324222875188 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, mq, beta, mu))| {
                assert_equal(
                    polarization_quark_l_thermal_part_landau_i(*q, *om, *p, *mq, *beta, *mu),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, mq, beta, mu))| {
                assert_equal(
                    oneloop__gluon__polarization_quark_l_thermal_part_landau_i(
                        *q, *om, *p, *mq, *beta, *mu,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_polarization_quark_t_thermal_part_landau_i() {
        use super::{
            ffi::oneloop__gluon__polarization_quark_t_thermal_part_landau_i,
            polarization_quark_t_thermal_part_landau_i,
        };

        assert_eq!(crate::consts::nc(), 3);

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 3.2, 0.8),
            (0.35, 0.21, 2.16, 1.2, 3.2, 0.8),
            (0.35, 0.75, 2.16, 1.2, 3.2, 0.8),
            (0.35, 0.75, 1.15, 1.2, 3.2, 0.8),
            (0.35, 0.75, 1.15, 3.76, 3.2, 0.8),
            (0.35, 0.75, 1.15, 3.76, 0.19, 0.8),
            (0.35, 0.75, 1.15, 3.76, 0.19, 5.9),
        ];

        let res: [C; 7] = [
            0.0001555678194717369 + 0. * I,
            1.801267697432198e-05 + 0. * I,
            -2.3742678845706716e-05 + 0. * I,
            -0.00015027646885955933 + 0. * I,
            -3.597043163726963e-08 + 0. * I,
            -0.0003218137257399576 + 0. * I,
            -0.00036053181120629027 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, mq, beta, mu))| {
                assert_equal(
                    polarization_quark_t_thermal_part_landau_i(*q, *om, *p, *mq, *beta, *mu),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, mq, beta, mu))| {
                assert_equal(
                    oneloop__gluon__polarization_quark_t_thermal_part_landau_i(
                        *q, *om, *p, *mq, *beta, *mu,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_zero_matsubara_polarization_glue_l_thermal_part_landau_i() {
        use super::zero_matsubara::{
            ffi::oneloop__zero_matsubara__gluon__polarization_glue_l_thermal_part_landau_i,
            polarization_glue_l_thermal_part_landau_i,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.0003208941584982947 + 0. * I,
            -0.0030075249405222253 + 0. * I,
            -0.01222669559611368 + 0. * I,
            0.0019273974485877518 + 0. * I,
            -18.665752022987885 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                polarization_glue_l_thermal_part_landau_i(*q, *p, *m, *beta),
                res[i],
            )
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__gluon__polarization_glue_l_thermal_part_landau_i(
                    *q, *p, *m, *beta,
                ),
                res[i],
            )
        });
    }

    #[test]
    fn test_zero_matsubara_polarization_glue_t_thermal_part_landau_i() {
        use super::zero_matsubara::{
            ffi::oneloop__zero_matsubara__gluon__polarization_glue_t_thermal_part_landau_i,
            polarization_glue_t_thermal_part_landau_i,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.007683815997018903 + 0. * I,
            -0.003429309709538143 + 0. * I,
            0.0029896620913042677 + 0. * I,
            -0.0006726588299993001 + 0. * I,
            1.046045679964218 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                polarization_glue_t_thermal_part_landau_i(*q, *p, *m, *beta),
                res[i],
            )
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__gluon__polarization_glue_t_thermal_part_landau_i(
                    *q, *p, *m, *beta,
                ),
                res[i],
            )
        });
    }

    #[test]
    fn test_zero_matsubara_polarization_quark_l_thermal_part_landau_i() {
        use super::zero_matsubara::{
            ffi::oneloop__zero_matsubara__gluon__polarization_quark_l_thermal_part_landau_i,
            polarization_quark_l_thermal_part_landau_i,
        };

        assert_eq!(crate::consts::nc(), 3);

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 3.2, 0.8),
            (0.35, 2.16, 1.2, 3.2, 0.8),
            (0.35, 1.15, 1.2, 3.2, 0.8),
            (0.35, 1.15, 3.76, 3.2, 0.8),
            (0.35, 1.15, 3.76, 0.19, 0.8),
            (0.35, 1.15, 3.76, 0.19, 5.9),
        ];

        let res: [R; 6] = [
            -0.0011644036243622508,
            -0.0004318279074659662,
            -0.001700976124595237,
            -2.0103148436168134e-06,
            -0.01798551978076218,
            -0.020149395452711206,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, mq, beta, mu))| {
                assert_equal(
                    polarization_quark_l_thermal_part_landau_i(*q, *p, *mq, *beta, *mu),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, mq, beta, mu))| {
                assert_equal(
                    oneloop__zero_matsubara__gluon__polarization_quark_l_thermal_part_landau_i(
                        *q, *p, *mq, *beta, *mu,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_zero_matsubara_polarization_quark_t_thermal_part_landau_i() {
        use super::zero_matsubara::{
            ffi::oneloop__zero_matsubara__gluon__polarization_quark_t_thermal_part_landau_i,
            polarization_quark_t_thermal_part_landau_i,
        };

        assert_eq!(crate::consts::nc(), 3);

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 3.2, 0.8),
            (0.35, 2.16, 1.2, 3.2, 0.8),
            (0.35, 1.15, 1.2, 3.2, 0.8),
            (0.35, 1.15, 3.76, 3.2, 0.8),
            (0.35, 1.15, 3.76, 0.19, 0.8),
            (0.35, 1.15, 3.76, 0.19, 5.9),
        ];

        let res: [R; 6] = [
            0.00018211145069647552,
            2.336631910538359e-05,
            9.444246149959117e-05,
            1.1909109425785625e-08,
            0.00010654625758190622,
            0.00011936506168260148,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, mq, beta, mu))| {
                assert_equal(
                    polarization_quark_t_thermal_part_landau_i(*q, *p, *mq, *beta, *mu),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, mq, beta, mu))| {
                assert_equal(
                    oneloop__zero_matsubara__gluon__polarization_quark_t_thermal_part_landau_i(
                        *q, *p, *mq, *beta, *mu,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_zero_momentum_polarization_glue_l_thermal_part_landau_i() {
        use super::zero_momentum::{
            ffi::oneloop__zero_momentum__gluon__polarization_glue_l_thermal_part_landau_i,
            polarization_glue_l_thermal_part_landau_i,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.003039308412097181 + 0. * I,
            -0.004004367287641101 + 0. * I,
            -0.0017951511301322243 + 0. * I,
            -0.0007127085951006912 + 0. * I,
            -0.28676848218619055 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                polarization_glue_l_thermal_part_landau_i(*q, *om, *m, *beta),
                res[i],
            )
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__gluon__polarization_glue_l_thermal_part_landau_i(
                    *q, *om, *m, *beta,
                ),
                res[i],
            )
        });
    }

    #[test]
    fn test_zero_momentum_polarization_glue_t_thermal_part_landau_i() {
        use super::zero_momentum::{
            ffi::oneloop__zero_momentum__gluon__polarization_glue_t_thermal_part_landau_i,
            polarization_glue_t_thermal_part_landau_i,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.003039308412097181 + 0. * I,
            -0.004004367287641101 + 0. * I,
            -0.0017951511301322243 + 0. * I,
            -0.0007127085951006912 + 0. * I,
            -0.28676848218619055 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                polarization_glue_t_thermal_part_landau_i(*q, *om, *m, *beta),
                res[i],
            )
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__gluon__polarization_glue_t_thermal_part_landau_i(
                    *q, *om, *m, *beta,
                ),
                res[i],
            )
        });
    }

    #[test]
    fn test_zero_momentum_polarization_quark_l_thermal_part_landau_i() {
        use super::zero_momentum::{
            ffi::oneloop__zero_momentum__gluon__polarization_quark_l_thermal_part_landau_i,
            polarization_quark_l_thermal_part_landau_i,
        };

        assert_eq!(crate::consts::nc(), 3);

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 3.2, 0.8),
            (0.35, 0.21, 1.2, 3.2, 0.8),
            (0.35, 0.75, 1.2, 3.2, 0.8),
            (0.35, 0.75, 3.76, 3.2, 0.8),
            (0.35, 0.75, 3.76, 0.19, 0.8),
            (0.35, 0.75, 3.76, 0.19, 5.9),
        ];

        let res: [R; 6] = [
            -0.0006552698065849495,
            -0.00030880656355125976,
            -0.0002853078006088784,
            -3.975991731807946e-08,
            -0.0003557168080792885,
            -0.00039851384461139094,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, mq, beta, mu))| {
                assert_equal(
                    polarization_quark_l_thermal_part_landau_i(*q, *om, *mq, *beta, *mu),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, mq, beta, mu))| {
                assert_equal(
                    oneloop__zero_momentum__gluon__polarization_quark_l_thermal_part_landau_i(
                        *q, *om, *mq, *beta, *mu,
                    ),
                    res[i],
                )
            });
    }

    #[test]
    fn test_zero_momentum_polarization_quark_t_thermal_part_landau_i() {
        use super::zero_momentum::{
            ffi::oneloop__zero_momentum__gluon__polarization_quark_t_thermal_part_landau_i,
            polarization_quark_t_thermal_part_landau_i,
        };

        assert_eq!(crate::consts::nc(), 3);

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 3.2, 0.8),
            (0.35, 0.21, 1.2, 3.2, 0.8),
            (0.35, 0.75, 1.2, 3.2, 0.8),
            (0.35, 0.75, 3.76, 3.2, 0.8),
            (0.35, 0.75, 3.76, 0.19, 0.8),
            (0.35, 0.75, 3.76, 0.19, 5.9),
        ];

        let res: [R; 6] = [
            -0.0006552698065849495,
            -0.00030880656355125976,
            -0.0002853078006088784,
            -3.975991731807946e-08,
            -0.0003557168080792885,
            -0.00039851384461139094,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, mq, beta, mu))| {
                assert_equal(
                    polarization_quark_t_thermal_part_landau_i(*q, *om, *mq, *beta, *mu),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, mq, beta, mu))| {
                assert_equal(
                    oneloop__zero_momentum__gluon__polarization_quark_t_thermal_part_landau_i(
                        *q, *om, *mq, *beta, *mu,
                    ),
                    res[i],
                )
            });
    }
}
