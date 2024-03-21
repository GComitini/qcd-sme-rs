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
    use crate::R;
    use peroxide::numerical::integral::Integral;

    pub fn thermal_polarization_l_landau_w_method(
        om: R,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> R {
        inlines::thermal_polarization_l_landau_w_method(om, p, m, beta, integral)
    }

    pub fn thermal_polarization_l_landau(om: R, p: R, m: R, beta: R) -> R {
        inlines::thermal_polarization_l_landau(om, p, m, beta)
    }

    pub fn thermal_polarization_t_landau_w_method(
        om: R,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> R {
        inlines::thermal_polarization_t_landau_w_method(om, p, m, beta, integral)
    }

    pub fn thermal_polarization_t_landau(om: R, p: R, m: R, beta: R) -> R {
        inlines::thermal_polarization_t_landau(om, p, m, beta)
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
pub(crate) mod ffi {}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::consts::get_default_integration_method;
    use crate::low_level::oneloop::thermal::*;
    use crate::{Num, C, R};
    use peroxide::numerical::integral::{integrate, Integral};

    #[inline(always)]
    pub fn thermal_polarization_l_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
        let s = om * om + p * p;
        let s2 = s * s;
        let s_inv = s.inv();
        let m2 = m * m;
        let m4 = m2 * m2;
        let a = s + m2 * 2.;
        let b = s + m2;

        let j_m_l_p_i = s_inv * (om * om * j_m_t_i(q, m, beta) + p * p * j_m_l_i(q, m, beta));
        let j_0_l_p_i = s_inv * (om * om * j_0_t_i(q, beta) + p * p * j_0_l_i(q, beta));
        let d_j_m_l_p_i = s_inv * (om * om * d_j_m_t_i(q, m, beta) + p * p * d_j_m_l_i(q, m, beta));
        let d2_j_m_l_p_i =
            s_inv * (om * om * d2_j_m_t_i(q, m, beta) + p * p * d2_j_m_l_i(q, m, beta));

        let t1 = (s2 * 3. / (m * 2.) - 1.) * i_0_0_l_i(q, om, p, beta);
        let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m2) + 4.)
            * i_m_m_l_i_same_mass(q, om, p, m, beta);
        let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_l_i(q, om, p, m, beta);
        let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, om, p, m, beta);
        let t5 = -b * s * 2. / m2 * i_m_0_i(q, om, p, m, beta);
        let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
        let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
        let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_l_i_same_mass(q, om, p, m, beta);
        let t9 = b * b / m2 * d_i_m_0_l_i(q, om, p, m, beta);
        let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, om, p, m, beta);
        let t11 = b * b * d_i_m_0_i(q, om, p, m, beta);
        let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

        let t13 = -m4 * d2_j_m_i(q, m, beta);
        let t14 = (j_m_l_p_i - j_0_l_p_i) / m2;
        let t15 = -d_j_m_l_p_i + d2_j_m_l_p_i * m2 / 2.;

        (t14 + t15 + t13) + t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11
    }

    #[inline(always)]
    pub fn thermal_polarization_t_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
        let s = om * om + p * p;
        let s2 = s * s;
        let m2 = m * m;
        let m4 = m2 * m2;
        let a = s + m2 * 2.;
        let b = s + m2;

        let t1 = (s2 * 3. / (m * 2.) - 1.) * i_0_0_t_i(q, om, p, beta);
        let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m2) + 4.)
            * i_m_m_t_i_same_mass(q, om, p, m, beta);
        let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_t_i(q, om, p, m, beta);
        let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, om, p, m, beta);
        let t5 = -b * s * 2. / m2 * i_m_0_i(q, om, p, m, beta);
        let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
        let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
        let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_t_i_same_mass(q, om, p, m, beta);
        let t9 = b * b / m2 * d_i_m_0_t_i(q, om, p, m, beta);
        let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, om, p, m, beta);
        let t11 = b * b * d_i_m_0_i(q, om, p, m, beta);
        let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

        let t13 = -m4 * d2_j_m_i(q, m, beta);
        let t14 = (j_m_t_i(q, m, beta) - j_0_t_i(q, beta)) / m2;
        let t15 = -d_j_m_t_i(q, m, beta) + d2_j_m_t_i(q, m, beta) * m2 / 2.;

        t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11 + t14 + t15 + t13
    }

    #[inline(always)]
    pub fn thermal_polarization_l_landau_w_method(
        om: R,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> R {
        integrate(
            |t| thermal_polarization_l_landau_i((1. - t) / t, om, p, m, beta).re / (t * t),
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn thermal_polarization_l_landau(om: R, p: R, m: R, beta: R) -> R {
        thermal_polarization_l_landau_w_method(om, p, m, beta, get_default_integration_method())
    }

    #[inline(always)]
    pub fn thermal_polarization_t_landau_w_method(
        om: R,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> R {
        integrate(
            |t| thermal_polarization_t_landau_i((1. - t) / t, om, p, m, beta).re / (t * t),
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn thermal_polarization_t_landau(om: R, p: R, m: R, beta: R) -> R {
        thermal_polarization_t_landau_w_method(om, p, m, beta, get_default_integration_method())
    }
}

pub mod zero_matsubara {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::{C, R};

        pub fn thermal_polarization_l_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::thermal_polarization_l_landau_i(q, p, m, beta)
        }

        pub fn thermal_polarization_t_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::thermal_polarization_t_landau_i(q, p, m, beta)
        }
    }

    mod ffi {}

    mod inlines {
        use crate::low_level::oneloop::thermal::zero_matsubara::*;
        use crate::low_level::oneloop::thermal::{
            d2_j_m_i, d2_j_m_l_i, d2_j_m_t_i, d_j_m_i, d_j_m_l_i, d_j_m_t_i, j_0_i, j_0_l_i,
            j_0_t_i, j_m_i, j_m_l_i, j_m_t_i,
        };
        use crate::{C, R};

        #[inline(always)]
        pub fn thermal_polarization_l_landau_i(q: R, p: R, m: R, beta: R) -> C {
            let s = p * p;
            let s2 = s * s;
            let m2 = m * m;
            let m4 = m2 * m2;
            let a = s + m2 * 2.;
            let b = s + m2;

            let j_m_l_p_i = j_m_l_i(q, m, beta);
            let j_0_l_p_i = j_0_l_i(q, beta);
            let d_j_m_l_p_i = d_j_m_l_i(q, m, beta);
            let d2_j_m_l_p_i = d2_j_m_l_i(q, m, beta);

            let t1 = (s2 * 3. / (m * 2.) - 1.) * i_0_0_l_i(q, p, beta);
            let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m2) + 4.)
                * i_m_m_l_i_same_mass(q, p, m, beta);
            let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_l_i(q, p, m, beta);
            let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, p, m, beta);
            let t5 = -b * s * 2. / m2 * i_m_0_i(q, p, m, beta);
            let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
            let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
            let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_l_i_same_mass(q, p, m, beta);
            let t9 = b * b / m2 * d_i_m_0_l_i(q, p, m, beta);
            let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, p, m, beta);
            let t11 = b * b * d_i_m_0_i(q, p, m, beta);
            let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

            let t13 = -m4 * d2_j_m_i(q, m, beta);
            let t14 = (j_m_l_p_i - j_0_l_p_i) / m2;
            let t15 = -d_j_m_l_p_i + d2_j_m_l_p_i * m2 / 2.;

            (t14 + t15 + t13) + t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11
        }

        #[inline(always)]
        pub fn thermal_polarization_t_landau_i(q: R, p: R, m: R, beta: R) -> C {
            let s = p * p;
            let s2 = s * s;
            let m2 = m * m;
            let m4 = m2 * m2;
            let a = s + m2 * 2.;
            let b = s + m2;

            let t1 = (s2 * 3. / (m * 2.) - 1.) * i_0_0_t_i(q, p, beta);
            let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m2) + 4.)
                * i_m_m_t_i_same_mass(q, p, m, beta);
            let t3 = -(s2 * 3. + s * m2 * 4. + m4) / m4 * i_m_0_t_i(q, p, m, beta);
            let t4 = a * s * 2. / m2 * i_m_m_i_same_mass(q, p, m, beta);
            let t5 = -b * s * 2. / m2 * i_m_0_i(q, p, m, beta);
            let t6 = -(s * 2. + m2 * 3.) / m2 * j_m_i(q, m, beta);
            let t7 = (s * 2. + m2) / m2 * j_0_i(q, beta);
            let t8 = -(a * a / m2 + m2 * 8.) * d_i_m_m_t_i_same_mass(q, p, m, beta);
            let t9 = b * b / m2 * d_i_m_0_t_i(q, p, m, beta);
            let t10 = -s * 2. * (s + m2 * 4.) * d_i_m_m_i_same_mass(q, p, m, beta);
            let t11 = b * b * d_i_m_0_i(q, p, m, beta);
            let t12 = (s + m2 * 3.) * d_j_m_i(q, m, beta);

            let t13 = -m4 * d2_j_m_i(q, m, beta);
            let t14 = (j_m_t_i(q, m, beta) - j_0_t_i(q, beta)) / m2;
            let t15 = -d_j_m_t_i(q, m, beta) + d2_j_m_t_i(q, m, beta) * m2 / 2.;

            t12 + (t7 + (t6 + t1)) + t2 + t3 + t4 + t5 + t8 + t9 + t10 + t11 + t14 + t15 + t13
        }
    }
}

pub mod zero_momentum {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::{Num, C, R};

        pub fn thermal_polarization_l_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            inlines::thermal_polarization_l_landau_i(q, om, m, beta)
        }

        pub fn thermal_polarization_t_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            inlines::thermal_polarization_t_landau_i(q, om, m, beta)
        }
    }

    mod ffi {}

    mod inlines {
        use crate::low_level::oneloop::thermal::zero_momentum::*;
        use crate::low_level::oneloop::thermal::{
            d2_j_m_i, d2_j_m_t_i, d_j_m_i, d_j_m_t_i, j_0_i, j_0_t_i, j_m_i, j_m_t_i,
        };
        use crate::{Num, C, R};

        #[inline(always)]
        pub fn thermal_polarization_l_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
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

            let t1 = (s2 * 3. / (m * 2.) - 1.) * i_0_0_l_i(q, om, beta);
            let t2 = ((s2 * 3. + s * m2 * 8. + 4. * m4) / (2. * m2) + 4.)
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
        pub fn thermal_polarization_t_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            thermal_polarization_l_landau_i(q, om, m, beta)
        }
    }
}
