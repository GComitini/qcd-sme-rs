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

    pub fn thermal_self_energy_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
        inlines::thermal_self_energy_landau_i(q, om, p, m, beta)
    }

    pub fn thermal_self_energy_landau_w_method(
        om: R,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> R {
        inlines::thermal_self_energy_landau_w_method(om, p, m, beta, integral)
    }

    pub fn thermal_self_energy_landau(om: R, p: R, m: R, beta: R) -> R {
        inlines::thermal_self_energy_landau(om, p, m, beta)
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
    pub fn thermal_self_energy_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
        let s = om * om + p * p;
        let m2 = m * m;
        let a = s + m2;

        let t1 = -a * s / (2. * m2) * i_m_0_i(q, om, p, m, beta);
        let t2 = s * s / (2. * m2) * i_0_0_i(q, om, p, beta);
        let t3 = (s * 2. - m2) / (4. * m2) * (j_m_i(q, m, beta) - j_0_i(q, beta));
        let t4 = a * a / 4. * d_i_m_0_i(q, om, p, m, beta);
        let t5 = -(s - m2) / 4. * d_j_m_i(q, m, beta);

        t5 + (t3 + t1) + t2 + t4
    }

    #[inline(always)]
    pub fn thermal_self_energy_landau_w_method(
        om: R,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> R {
        integrate(
            |t| thermal_self_energy_landau_i((1. - t) / t, om, p, m, beta).re / (t * t),
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn thermal_self_energy_landau(om: R, p: R, m: R, beta: R) -> R {
        thermal_self_energy_landau_w_method(om, p, m, beta, get_default_integration_method())
    }
}

pub mod zero_matsubara {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::Integral;
        use crate::{C, R};

        pub fn thermal_self_energy_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::thermal_self_energy_landau_i(q, p, m, beta)
        }

        pub fn thermal_self_energy_landau_w_method(p: R, m: R, beta: R, integral: Integral) -> R {
            inlines::thermal_self_energy_landau_w_method(p, m, beta, integral)
        }

        pub fn thermal_self_energy_landau(p: R, m: R, beta: R) -> R {
            inlines::thermal_self_energy_landau(p, m, beta)
        }
    }

    mod ffi {}

    mod inlines {
        use crate::consts::get_default_integration_method;
        use crate::low_level::oneloop::thermal::zero_matsubara::*;
        use crate::low_level::oneloop::thermal::{d_j_m_i, j_0_i, j_m_i};
        use crate::{C, R};
        use peroxide::numerical::integral::{integrate, Integral};

        #[inline(always)]
        pub fn thermal_self_energy_landau_i(q: R, p: R, m: R, beta: R) -> C {
            let s = p * p;
            let m2 = m * m;
            let a = s + m2;

            let t1 = -a * s / (2. * m2) * i_m_0_i(q, p, m, beta);
            let t2 = s * s / (2. * m2) * i_0_0_i(q, p, beta);
            let t3 = (s * 2. - m2) / (4. * m2) * (j_m_i(q, m, beta) - j_0_i(q, beta));
            let t4 = a * a / 4. * d_i_m_0_i(q, p, m, beta);
            let t5 = -(s - m2) / 4. * d_j_m_i(q, m, beta);

            t5 + (t3 + t1) + t2 + t4
        }

        #[inline(always)]
        pub fn thermal_self_energy_landau_w_method(p: R, m: R, beta: R, integral: Integral) -> R {
            integrate(
                |t| thermal_self_energy_landau_i((1. - t) / t, p, m, beta).re / (t * t),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn thermal_self_energy_landau(p: R, m: R, beta: R) -> R {
            thermal_self_energy_landau_w_method(p, m, beta, get_default_integration_method())
        }
    }
}

pub mod zero_momentum {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::Integral;
        use crate::{C, R};

        pub fn thermal_self_energy_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::thermal_self_energy_landau_i(q, p, m, beta)
        }

        pub fn thermal_self_energy_landau_w_method(om: R, m: R, beta: R, integral: Integral) -> R {
            inlines::thermal_self_energy_landau_w_method(om, m, beta, integral)
        }

        pub fn thermal_self_energy_landau(om: R, m: R, beta: R) -> R {
            inlines::thermal_self_energy_landau(om, m, beta)
        }
    }

    mod ffi {}

    mod inlines {
        use crate::consts::get_default_integration_method;
        use crate::low_level::oneloop::thermal::zero_momentum::*;
        use crate::low_level::oneloop::thermal::{d_j_m_i, j_0_i, j_m_i};
        use crate::{Num, C, R};
        use peroxide::numerical::integral::{integrate, Integral};

        #[inline(always)]
        pub fn thermal_self_energy_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
            let s = om * om;
            let m2 = m * m;
            let a = s + m2;

            let t1 = -a * s / (2. * m2) * i_m_0_i(q, om, m, beta);
            let t2 = s * s / (2. * m2) * i_0_0_i(q, om, beta);
            let t3 = (s * 2. - m2) / (4. * m2) * (j_m_i(q, m, beta) - j_0_i(q, beta));
            let t4 = a * a / 4. * d_i_m_0_i(q, om, m, beta);
            let t5 = -(s - m2) / 4. * d_j_m_i(q, m, beta);

            t5 + (t3 + t1) + t2 + t4
        }

        #[inline(always)]
        pub fn thermal_self_energy_landau_w_method(om: R, m: R, beta: R, integral: Integral) -> R {
            integrate(
                |t| thermal_self_energy_landau_i((1. - t) / t, om, m, beta).re / (t * t),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn thermal_self_energy_landau(om: R, m: R, beta: R) -> R {
            thermal_self_energy_landau_w_method(om, m, beta, get_default_integration_method())
        }
    }
}
