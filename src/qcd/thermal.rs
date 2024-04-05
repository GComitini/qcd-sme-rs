pub mod ghost {
    pub use crate::ym::thermal::ghost::*;
}

pub mod gluon {
    use crate::{Num, C, R};

    pub fn dressing_l_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
        inlines::dressing_l_inv_landau(om, p, m, beta, mu, f0)
    }

    pub fn dressing_t_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
        inlines::dressing_t_inv_landau(om, p, m, beta, mu, f0)
    }

    pub fn dressing_l_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
        inlines::dressing_l_landau(om, p, m, beta, mu, f0)
    }

    pub fn dressing_t_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
        inlines::dressing_t_landau(om, p, m, beta, mu, f0)
    }

    pub fn propagator_l_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
        inlines::propagator_l_landau(om, p, m, beta, mu, f0)
    }

    pub fn propagator_t_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
        inlines::propagator_t_landau(om, p, m, beta, mu, f0)
    }

    pub(crate) mod ffi {}

    pub(crate) mod inlines {
        use crate::consts::{m_quark, nf};
        use crate::low_level::oneloop::thermal::gluon as gluon_thermal_parts;
        use crate::qcd::gluon as gluon_vac;
        use crate::{Num, C, R};
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        #[inline(always)]
        pub fn dressing_l_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, f0)
                - sdim.inv()
                    * PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(om, p, m, beta)
                        + (nf() as R)
                            * gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                                om,
                                p,
                                m_quark(),
                                beta,
                                mu,
                            ))
        }

        #[inline(always)]
        pub fn dressing_t_inv_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
            let sdim = om * om + p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, f0)
                - sdim.inv()
                    * PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(om, p, m, beta)
                        + (nf() as R)
                            * gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                                om,
                                p,
                                m_quark(),
                                beta,
                                mu,
                            ))
        }

        #[inline(always)]
        pub fn dressing_l_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
            dressing_l_inv_landau(om, p, m, beta, mu, f0).inv()
        }

        #[inline(always)]
        pub fn dressing_t_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
            dressing_t_inv_landau(om, p, m, beta, mu, f0).inv()
        }

        #[inline(always)]
        pub fn propagator_l_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
            (om * om + p * p).inv() * dressing_l_landau(om, p, m, beta, mu, f0)
        }

        #[inline(always)]
        pub fn propagator_t_landau<T: Num>(om: T, p: R, m: R, beta: R, mu: R, f0: R) -> C {
            (om * om + p * p).inv() * dressing_t_landau(om, p, m, beta, mu, f0)
        }
    }
}

pub mod zero_matsubara {
    pub mod ghost {
        pub use crate::ym::thermal::zero_matsubara::ghost::*;
    }

    pub mod gluon {
        use crate::R;

        pub fn dressing_l_inv_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
            inlines::dressing_l_inv_landau(p, m, beta, mu, f0)
        }

        pub fn dressing_t_inv_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
            inlines::dressing_t_inv_landau(p, m, beta, mu, f0)
        }

        pub fn dressing_l_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
            inlines::dressing_l_landau(p, m, beta, mu, f0)
        }

        pub fn dressing_t_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
            inlines::dressing_t_landau(p, m, beta, mu, f0)
        }

        pub fn propagator_l_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
            inlines::propagator_l_landau(p, m, beta, mu, f0)
        }

        pub fn propagator_t_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
            inlines::propagator_t_landau(p, m, beta, mu, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::consts::{m_quark, nf};
            use crate::low_level::oneloop::thermal::gluon::zero_matsubara as gluon_thermal_parts;
            use crate::qcd::gluon as gluon_vac;
            use crate::R;
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_l_inv_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - PREFACTOR
                        * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(p, m, beta)
                            + (nf() as R)
                                * gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                                    p,
                                    m_quark(),
                                    beta,
                                    mu,
                                ))
                        / sdim
            }

            #[inline(always)]
            pub fn dressing_t_inv_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
                let sdim = p * p;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - PREFACTOR
                        * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(p, m, beta)
                            + (nf() as R)
                                * gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                                    p,
                                    m_quark(),
                                    beta,
                                    mu,
                                ))
                        / sdim
            }

            #[inline(always)]
            pub fn dressing_l_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
                1. / dressing_l_inv_landau(p, m, beta, mu, f0)
            }

            #[inline(always)]
            pub fn dressing_t_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
                1. / dressing_t_inv_landau(p, m, beta, mu, f0)
            }

            #[inline(always)]
            pub fn propagator_l_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
                dressing_l_landau(p, m, beta, mu, f0) / (p * p)
            }

            #[inline(always)]
            pub fn propagator_t_landau(p: R, m: R, beta: R, mu: R, f0: R) -> R {
                dressing_t_landau(p, m, beta, mu, f0) / (p * p)
            }
        }
    }
}

pub mod zero_momentum {
    pub mod ghost {
        pub use crate::ym::thermal::zero_momentum::ghost::*;
    }

    pub mod gluon {
        use crate::{Num, C, R};

        pub fn dressing_l_inv_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
            inlines::dressing_l_inv_landau(om, m, beta, mu, f0)
        }

        pub fn dressing_t_inv_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
            inlines::dressing_t_inv_landau(om, m, beta, mu, f0)
        }

        pub fn dressing_l_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
            inlines::dressing_l_landau(om, m, beta, mu, f0)
        }

        pub fn dressing_t_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
            inlines::dressing_t_landau(om, m, beta, mu, f0)
        }

        pub fn propagator_l_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
            inlines::propagator_l_landau(om, m, beta, mu, f0)
        }

        pub fn propagator_t_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
            inlines::propagator_t_landau(om, m, beta, mu, f0)
        }

        pub(crate) mod ffi {}

        pub(crate) mod inlines {
            use crate::consts::{m_quark, nf};
            use crate::low_level::oneloop::thermal::gluon::zero_momentum as gluon_thermal_parts;
            use crate::qcd::ghost as gluon_vac;
            use crate::{Num, C, R};
            use std::f64::consts::PI;

            const PREFACTOR: R = (16. * PI * PI) / 3.;

            #[inline(always)]
            pub fn dressing_l_inv_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
                let sdim = om * om;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - sdim.inv()
                        * PREFACTOR
                        * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(
                            om, m, beta,
                        ) + (nf() as R)
                            * gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                                om,
                                m_quark(),
                                beta,
                                mu,
                            ))
            }

            #[inline(always)]
            pub fn dressing_t_inv_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
                let sdim = om * om;
                let s = sdim / (m * m);
                gluon_vac::dressing_inv_landau(s, f0)
                    - sdim.inv()
                        * PREFACTOR
                        * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(
                            om, m, beta,
                        ) + (nf() as R)
                            * gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                                om,
                                m_quark(),
                                beta,
                                mu,
                            ))
            }

            #[inline(always)]
            pub fn dressing_l_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
                dressing_l_inv_landau(om, m, beta, mu, f0).inv()
            }

            #[inline(always)]
            pub fn dressing_t_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
                dressing_t_inv_landau(om, m, beta, mu, f0).inv()
            }

            #[inline(always)]
            pub fn propagator_l_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
                (om * om).inv() * dressing_l_landau(om, m, beta, mu, f0)
            }

            #[inline(always)]
            pub fn propagator_t_landau<T: Num>(om: T, m: R, beta: R, mu: R, f0: R) -> C {
                (om * om).inv() * dressing_t_landau(om, m, beta, mu, f0)
            }
        }
    }
}
