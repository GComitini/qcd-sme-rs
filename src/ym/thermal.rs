pub mod ghost {
    use crate::R;

    pub fn dressing_inv_landau(om: R, p: R, m: R, beta: R, g0: R) -> R {
        inlines::dressing_inv_landau(om, p, m, beta, g0)
    }

    pub(crate) mod ffi {}

    pub(crate) mod inlines {
        use crate::low_level::oneloop::thermal::ghost::self_energy_thermal_part_landau;
        use crate::ym::ghost as ghost_vac;
        use crate::R;
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        #[inline(always)]
        pub fn dressing_inv_landau(om: R, p: R, m: R, beta: R, g0: R) -> R {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            ghost_vac::dressing_inv_landau(s, g0)
                + PREFACTOR * self_energy_thermal_part_landau(om, p, m, beta) / sdim
        }
    }
}

pub mod gluon {
    use crate::R;

    pub fn dressing_l_inv_landau(om: R, p: R, m: R, beta: R, f0: R) -> R {
        inlines::dressing_l_inv_landau(om, p, m, beta, f0)
    }

    pub fn dressing_t_inv_landau(om: R, p: R, m: R, beta: R, f0: R) -> R {
        inlines::dressing_t_inv_landau(om, p, m, beta, f0)
    }

    pub(crate) mod ffi {}

    pub(crate) mod inlines {
        use crate::low_level::oneloop::thermal::gluon::{
            polarization_l_thermal_part_landau, polarization_t_thermal_part_landau,
        };
        use crate::ym::gluon as gluon_vac;
        use crate::R;
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        #[inline(always)]
        pub fn dressing_l_inv_landau(om: R, p: R, m: R, beta: R, f0: R) -> R {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, f0)
                - PREFACTOR * polarization_l_thermal_part_landau(om, p, m, beta) / sdim
        }

        #[inline(always)]
        pub fn dressing_t_inv_landau(om: R, p: R, m: R, beta: R, f0: R) -> R {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, f0)
                - PREFACTOR * polarization_t_thermal_part_landau(om, p, m, beta) / sdim
        }
    }
}

pub mod zero_matsubara {
    pub mod ghost {
        use crate::R;

        pub fn dressing_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_inv_landau(p, m, beta, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::ghost::zero_matsubara::self_energy_thermal_part_landau;
            use crate::ym::ghost as ghost_vac;
            use crate::R;
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                ghost_vac::dressing_inv_landau(s, f0)
                    + PREFACTOR * self_energy_thermal_part_landau(p, m, beta) / sdim
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

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::gluon::zero_matsubara::{
                polarization_l_thermal_part_landau, polarization_t_thermal_part_landau,
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
                    - PREFACTOR * polarization_l_thermal_part_landau(p, m, beta) / sdim
            }

            #[inline(always)]
            pub fn dressing_t_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - PREFACTOR * polarization_t_thermal_part_landau(p, m, beta) / sdim
            }
        }
    }
}

pub mod zero_momentum {
    pub mod ghost {
        use crate::R;

        pub fn dressing_inv_landau(p: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_inv_landau(p, m, beta, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::ghost::zero_momentum::self_energy_thermal_part_landau;
            use crate::ym::ghost as ghost_vac;
            use crate::R;
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_inv_landau(om: R, m: R, beta: R, f0: R) -> R {
                let sdim = om * om;
                let s = sdim / (m * m);
                ghost_vac::dressing_inv_landau(s, f0)
                    + PREFACTOR * self_energy_thermal_part_landau(om, m, beta) / sdim
            }
        }
    }

    pub mod gluon {
        use crate::R;

        pub fn dressing_l_inv_landau(om: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_l_inv_landau(om, m, beta, f0)
        }

        pub fn dressing_t_inv_landau(om: R, m: R, beta: R, f0: R) -> R {
            inlines::dressing_t_inv_landau(om, m, beta, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::low_level::oneloop::thermal::gluon::zero_momentum::polarization_l_thermal_part_landau;
            use crate::ym::gluon as gluon_vac;
            use crate::R;
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_l_inv_landau(om: R, m: R, beta: R, f0: R) -> R {
                let sdim = om * om;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - PREFACTOR * polarization_l_thermal_part_landau(om, m, beta) / sdim
            }

            #[inline(always)]
            pub fn dressing_t_inv_landau(om: R, m: R, beta: R, f0: R) -> R {
                dressing_l_inv_landau(om, m, beta, f0)
            }
        }
    }
}
