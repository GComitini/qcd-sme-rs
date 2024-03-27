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

    pub fn self_energy_thermal_part_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
        inlines::self_energy_thermal_part_landau_i(q, om, p, m, beta)
    }

    pub fn self_energy_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> C {
        inlines::self_energy_thermal_part_landau_w_method(om, p, m, beta, integral)
    }

    pub fn self_energy_thermal_part_landau<T: Num>(om: T, p: R, m: R, beta: R) -> C {
        inlines::self_energy_thermal_part_landau(om, p, m, beta)
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
    pub extern "C" fn oneloop__ghost__self_energy_thermal_part_landau_i(
        q: R,
        om: R,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::self_energy_thermal_part_landau_i(q, om, p, m, beta)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::consts::get_default_integration_method;
    use crate::low_level::oneloop::thermal::*;
    use crate::{Num, C, R};
    use peroxide::numerical::integral::{complex_integrate, Integral};

    #[inline(always)]
    pub fn self_energy_thermal_part_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
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
    pub fn self_energy_thermal_part_landau_w_method<T: Num>(
        om: T,
        p: R,
        m: R,
        beta: R,
        integral: Integral,
    ) -> C {
        complex_integrate(
            |t| self_energy_thermal_part_landau_i((1. - t) / t, om, p, m, beta) / (t * t),
            (0., 1.),
            integral,
        )
    }

    #[inline(always)]
    pub fn self_energy_thermal_part_landau<T: Num>(om: T, p: R, m: R, beta: R) -> C {
        self_energy_thermal_part_landau_w_method(om, p, m, beta, get_default_integration_method())
    }
}

pub mod zero_matsubara {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::Integral;
        use crate::{C, R};

        pub fn self_energy_thermal_part_landau_i(q: R, p: R, m: R, beta: R) -> C {
            inlines::self_energy_thermal_part_landau_i(q, p, m, beta)
        }

        pub fn self_energy_thermal_part_landau_w_method(
            p: R,
            m: R,
            beta: R,
            integral: Integral,
        ) -> R {
            inlines::self_energy_thermal_part_landau_w_method(p, m, beta, integral)
        }

        pub fn self_energy_thermal_part_landau(p: R, m: R, beta: R) -> R {
            inlines::self_energy_thermal_part_landau(p, m, beta)
        }
    }

    pub(crate) mod ffi {
        use super::inlines;
        use crate::{C, R};

        #[no_mangle]
        pub extern "C" fn oneloop__zero_matsubara__ghost__self_energy_thermal_part_landau_i(
            q: R,
            p: R,
            m: R,
            beta: R,
        ) -> C {
            inlines::self_energy_thermal_part_landau_i(q, p, m, beta)
        }
    }

    mod inlines {
        use crate::consts::get_default_integration_method;
        use crate::low_level::oneloop::thermal::zero_matsubara::*;
        use crate::low_level::oneloop::thermal::{d_j_m_i, j_0_i, j_m_i};
        use crate::{C, R};
        use peroxide::numerical::integral::{integrate, Integral};

        #[inline(always)]
        pub fn self_energy_thermal_part_landau_i(q: R, p: R, m: R, beta: R) -> C {
            let s = p * p;
            let m2 = m * m;
            let a = s + m2;

            let t1 = -a * s / (2. * m2) * i_m_0_i(q, p, m, beta);
            let t2 = s * s / (2. * m2) * i_0_0_i(q, p, beta);
            let t3 = s.mul_add(2., -m2) / (4. * m2) * (j_m_i(q, m, beta) - j_0_i(q, beta));
            let t4 = a * a / 4. * d_i_m_0_i(q, p, m, beta);
            let t5 = -(s - m2) / 4. * d_j_m_i(q, m, beta);

            t5 + (t3 + t1) + t2 + t4
        }

        #[inline(always)]
        pub fn self_energy_thermal_part_landau_w_method(
            p: R,
            m: R,
            beta: R,
            integral: Integral,
        ) -> R {
            integrate(
                |t| self_energy_thermal_part_landau_i((1. - t) / t, p, m, beta).re / (t * t),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn self_energy_thermal_part_landau(p: R, m: R, beta: R) -> R {
            self_energy_thermal_part_landau_w_method(p, m, beta, get_default_integration_method())
        }
    }
}

pub mod zero_momentum {
    pub use native::*;

    mod native {
        use super::inlines;
        use crate::Integral;
        use crate::{Num, C, R};

        pub fn self_energy_thermal_part_landau_i(q: R, om: R, m: R, beta: R) -> C {
            inlines::self_energy_thermal_part_landau_i(q, om, m, beta)
        }

        pub fn self_energy_thermal_part_landau_w_method<T: Num>(
            om: T,
            m: R,
            beta: R,
            integral: Integral,
        ) -> C {
            inlines::self_energy_thermal_part_landau_w_method(om, m, beta, integral)
        }

        pub fn self_energy_thermal_part_landau<T: Num>(om: T, m: R, beta: R) -> C {
            inlines::self_energy_thermal_part_landau(om, m, beta)
        }
    }

    pub(crate) mod ffi {
        use super::inlines;
        use crate::{C, R};

        #[no_mangle]
        pub extern "C" fn oneloop__zero_momentum__ghost__self_energy_thermal_part_landau_i(
            q: R,
            om: R,
            m: R,
            beta: R,
        ) -> C {
            inlines::self_energy_thermal_part_landau_i(q, om, m, beta)
        }
    }

    mod inlines {
        use crate::consts::get_default_integration_method;
        use crate::low_level::oneloop::thermal::zero_momentum::*;
        use crate::low_level::oneloop::thermal::{d_j_m_i, j_0_i, j_m_i};
        use crate::{Num, C, R};
        use peroxide::numerical::integral::{complex_integrate, Integral};

        #[inline(always)]
        pub fn self_energy_thermal_part_landau_i<T: Num>(q: R, om: T, m: R, beta: R) -> C {
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
        pub fn self_energy_thermal_part_landau_w_method<T: Num>(
            om: T,
            m: R,
            beta: R,
            integral: Integral,
        ) -> C {
            complex_integrate(
                |t| self_energy_thermal_part_landau_i((1. - t) / t, om, m, beta) / (t * t),
                (0., 1.),
                integral,
            )
        }

        #[inline(always)]
        pub fn self_energy_thermal_part_landau<T: Num>(om: T, m: R, beta: R) -> C {
            self_energy_thermal_part_landau_w_method(om, m, beta, get_default_integration_method())
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
    fn test_self_energy_thermal_part_landau_i() {
        use super::{
            ffi::oneloop__ghost__self_energy_thermal_part_landau_i,
            self_energy_thermal_part_landau_i,
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
            1.0169819975930068e-05 + 0. * I,
            -0.00049390422718809 + 0. * I,
            -0.0002989874313323077 + 2.9161839440024374e-20 * I,
            -0.00036085673585223516 + 0. * I,
            -3.7033632854858595e-05 + 0. * I,
            0.0071054002993155285 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    self_energy_thermal_part_landau_i(*q, *om, *p, *m, *beta),
                    res[i],
                )
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    oneloop__ghost__self_energy_thermal_part_landau_i(*q, *om, *p, *m, *beta),
                    res[i],
                )
            });
    }

    #[test]
    fn test_zero_matsubara_self_energy_thermal_part_landau_i() {
        use super::zero_matsubara::{
            ffi::oneloop__zero_matsubara__ghost__self_energy_thermal_part_landau_i,
            self_energy_thermal_part_landau_i,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.0008113805733283437 + 0. * I,
            -0.0005651397812783704 + 0. * I,
            0.0005722445733281894 + 0. * I,
            0.00012338612611931306 + 0. * I,
            0.013490773390605788 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(self_energy_thermal_part_landau_i(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__ghost__self_energy_thermal_part_landau_i(
                    *q, *p, *m, *beta,
                ),
                res[i],
            )
        });
    }

    #[test]
    fn test_zero_momentum_self_energy_thermal_part_landau_i() {
        use super::zero_momentum::{
            ffi::oneloop__zero_momentum__ghost__self_energy_thermal_part_landau_i,
            self_energy_thermal_part_landau_i,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -4.858855776705466e-05 + 0. * I,
            -0.00018227434215185073 + 0. * I,
            -0.0010742009005355816 + 0. * I,
            -0.00015064888460959 + 0. * I,
            -0.004608109313031232 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                self_energy_thermal_part_landau_i(*q, *om, *m, *beta),
                res[i],
            )
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__ghost__self_energy_thermal_part_landau_i(
                    *q, *om, *m, *beta,
                ),
                res[i],
            )
        });
    }
}
