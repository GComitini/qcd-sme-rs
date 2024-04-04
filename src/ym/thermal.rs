pub mod ghost {
    use crate::{C, R};

    pub fn dressing_inv_landau(om: R, p: R, m: R, beta: R, g0: R) -> C {
        inlines::dressing_inv_landau(om, p, m, beta, g0)
    }

    pub fn dressing_landau(om: R, p: R, m: R, beta: R, g0: R) -> C {
        inlines::dressing_landau(om, p, m, beta, g0)
    }

    pub fn propagator_landau(om: R, p: R, m: R, beta: R, g0: R) -> C {
        inlines::propagator_landau(om, p, m, beta, g0)
    }

    pub(crate) mod ffi {}

    pub(crate) mod inlines {
        use crate::low_level::oneloop::thermal::ghost::self_energy_thermal_part_landau;
        use crate::ym::ghost as ghost_vac;
        use crate::{Num, C, R};
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        #[inline(always)]
        pub fn dressing_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, g0: R) -> C {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            ghost_vac::dressing_inv_landau(s, g0)
                + sdim.inv() * PREFACTOR * self_energy_thermal_part_landau(om, p, m, beta)
        }

        #[inline(always)]
        pub fn dressing_landau<T: Num>(om: T, p: R, m: R, beta: R, g0: R) -> C {
            dressing_inv_landau(om, p, m, beta, g0).inv()
        }

        #[inline(always)]
        pub fn propagator_landau<T: Num>(om: T, p: R, m: R, beta: R, g0: R) -> C {
            (om * om + p * p).inv() * dressing_landau(om, p, m, beta, g0)
        }
    }
}

pub mod gluon {
    use crate::{Num, C, R};

    pub fn dressing_l_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
        inlines::dressing_l_inv_landau(om, p, m, beta, f0)
    }

    pub fn dressing_t_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
        inlines::dressing_t_inv_landau(om, p, m, beta, f0)
    }

    pub fn dressing_l_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
        inlines::dressing_l_landau(om, p, m, beta, f0)
    }

    pub fn dressing_t_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
        inlines::dressing_t_landau(om, p, m, beta, f0)
    }

    pub fn propagator_l_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
        inlines::propagator_l_landau(om, p, m, beta, f0)
    }

    pub fn propagator_t_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
        inlines::propagator_t_landau(om, p, m, beta, f0)
    }

    pub(crate) mod ffi {}

    pub(crate) mod inlines {
        use crate::low_level::oneloop::thermal::gluon::{
            polarization_glue_l_thermal_part_landau, polarization_glue_t_thermal_part_landau,
        };
        use crate::ym::gluon as gluon_vac;
        use crate::{Num, C, R};
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        #[inline(always)]
        pub fn dressing_l_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, f0)
                - sdim.inv() * PREFACTOR * polarization_glue_l_thermal_part_landau(om, p, m, beta)
        }

        #[inline(always)]
        pub fn dressing_t_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, f0)
                - sdim.inv() * PREFACTOR * polarization_glue_t_thermal_part_landau(om, p, m, beta)
        }

        #[inline(always)]
        pub fn dressing_l_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
            dressing_l_inv_landau(om, p, m, beta, f0).inv()
        }

        #[inline(always)]
        pub fn dressing_t_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
            dressing_t_inv_landau(om, p, m, beta, f0).inv()
        }

        #[inline(always)]
        pub fn propagator_l_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
            (om * om + p * p).inv() * dressing_l_landau(om, p, m, beta, f0)
        }

        #[inline(always)]
        pub fn propagator_t_landau<T: Num>(om: T, p: R, m: R, beta: R, f0: R) -> C {
            (om * om + p * p).inv() * dressing_t_landau(om, p, m, beta, f0)
        }
    }
}

pub mod zero_matsubara {
    pub mod ghost {
        use crate::R;

        pub fn dressing_inv_landau(p: R, m: R, beta: R, g0: R) -> R {
            inlines::dressing_inv_landau(p, m, beta, g0)
        }

        pub fn dressing_landau(p: R, m: R, beta: R, g0: R) -> R {
            inlines::dressing_landau(p, m, beta, g0)
        }

        pub fn propagator_landau(p: R, m: R, beta: R, g0: R) -> R {
            inlines::propagator_landau(p, m, beta, g0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::ghost::zero_matsubara::self_energy_thermal_part_landau;
            use crate::ym::ghost as ghost_vac;
            use crate::R;
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_inv_landau(p: R, m: R, beta: R, g0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                ghost_vac::dressing_inv_landau(s, g0)
                    + PREFACTOR * self_energy_thermal_part_landau(p, m, beta) / sdim
            }

            #[inline(always)]
            pub fn dressing_landau(p: R, m: R, beta: R, g0: R) -> R {
                1. / dressing_inv_landau(p, m, beta, g0)
            }

            #[inline(always)]
            pub fn propagator_landau(p: R, m: R, beta: R, g0: R) -> R {
                dressing_landau(p, m, beta, g0) / (p * p)
            }
        }
    }

    pub mod gluon {
        use crate::R;

        pub fn dressing_l_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_l_inv_landau(p, m, beta, f0)
        }

        pub fn dressing_t_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_t_inv_landau(p, m, beta, f0)
        }

        pub fn dressing_l_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_l_landau(p, m, beta, f0)
        }

        pub fn dressing_t_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_t_landau(p, m, beta, f0)
        }

        pub fn propagator_l_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::propagator_l_landau(p, m, beta, f0)
        }

        pub fn propagator_t_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::propagator_t_landau(p, m, beta, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::gluon::zero_matsubara::{
                polarization_glue_l_thermal_part_landau, polarization_glue_t_thermal_part_landau,
            };
            use crate::ym::gluon as gluon_vac;
            use crate::R;
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_l_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - PREFACTOR * polarization_glue_l_thermal_part_landau(p, m, beta) / sdim
            }

            #[inline(always)]
            pub fn dressing_t_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - PREFACTOR * polarization_glue_t_thermal_part_landau(p, m, beta) / sdim
            }

            #[inline(always)]
            pub fn dressing_l_landau(p: R, m: R, beta: R, f0: R) -> R {
                1. / dressing_l_inv_landau(p, m, beta, f0)
            }

            #[inline(always)]
            pub fn dressing_t_landau(p: R, m: R, beta: R, f0: R) -> R {
                1. / dressing_t_inv_landau(p, m, beta, f0)
            }

            #[inline(always)]
            pub fn propagator_l_landau(p: R, m: R, beta: R, f0: R) -> R {
                dressing_l_landau(p, m, beta, f0) / (p * p)
            }

            #[inline(always)]
            pub fn propagator_t_landau(p: R, m: R, beta: R, f0: R) -> R {
                dressing_t_landau(p, m, beta, f0) / (p * p)
            }
        }
    }
}

pub mod zero_momentum {
    pub mod ghost {
        use crate::{Num, C, R};

        pub fn dressing_inv_landau<T: Num>(om: R, m: R, beta: R, g0: R) -> C {
            inlines::dressing_inv_landau(om, m, beta, g0)
        }

        pub fn dressing_landau<T: Num>(om: R, m: R, beta: R, g0: R) -> C {
            inlines::dressing_landau(om, m, beta, g0)
        }

        pub fn propagator_landau<T: Num>(om: R, m: R, beta: R, g0: R) -> C {
            inlines::propagator_landau(om, m, beta, g0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::ghost::zero_momentum::self_energy_thermal_part_landau;
            use crate::ym::ghost as ghost_vac;
            use crate::{Num, C, R};
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_inv_landau<T: Num>(om: T, m: R, beta: R, g0: R) -> C {
                let sdim = om * om;
                let s = sdim / (m * m);
                ghost_vac::dressing_inv_landau(s, g0)
                    + sdim.inv() * PREFACTOR * self_energy_thermal_part_landau(om, m, beta)
            }

            #[inline(always)]
            pub fn dressing_landau<T: Num>(om: T, m: R, beta: R, g0: R) -> C {
                1. / dressing_inv_landau(om, m, beta, g0)
            }

            #[inline(always)]
            pub fn propagator_landau<T: Num>(om: T, m: R, beta: R, g0: R) -> C {
                (om * om).inv() * dressing_landau(om, m, beta, g0)
            }
        }
    }

    pub mod gluon {
        use crate::{Num, C, R};

        pub fn dressing_l_inv_landau<T: Num>(om: T, m: R, beta: R, f0: R) -> C {
            inlines::dressing_l_inv_landau(om, m, beta, f0)
        }

        pub fn dressing_t_inv_landau<T: Num>(om: T, m: R, beta: R, f0: R) -> C {
            inlines::dressing_t_inv_landau(om, m, beta, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::gluon::zero_momentum::polarization_glue_l_thermal_part_landau;
            use crate::ym::gluon as gluon_vac;
            use crate::{Num, C, R};
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_l_inv_landau<T: Num>(om: T, m: R, beta: R, f0: R) -> C {
                let sdim = om * om;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - sdim.inv() * PREFACTOR * polarization_glue_l_thermal_part_landau(om, m, beta)
            }

            #[inline(always)]
            pub fn dressing_t_inv_landau<T: Num>(om: T, m: R, beta: R, f0: R) -> C {
                dressing_l_inv_landau(om, m, beta, f0)
            }
        }
    }
}
