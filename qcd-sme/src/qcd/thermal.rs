// At variance with the rest of crate::qcd, in this module quark masses are dimensionful.

pub mod ghost {
    pub use crate::ym::thermal::ghost::*;
}

pub mod gluon {
    use crate::consts::nf;
    use crate::low_level::oneloop::thermal::gluon as gluon_thermal_parts;
    use crate::qcd::gluon as gluon_vac;
    use crate::qcd::FieldConfig;
    use crate::{Num, C, R};
    use std::f64::consts::PI;

    const PREFACTOR: R = (16. * PI * PI) / 3.;

    pub fn dressing_l_inv_landau<T: Num>(om: T, p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
        let sdim = om * om + p * p;
        let s = sdim / (m * m);
        gluon_vac::dressing_inv_landau(s, mq / m, f0)
            - sdim.inv()
                * PREFACTOR
                * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(om, p, m, beta)
                    + (nf() as R)
                        * gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                            om, p, mq, beta, mu,
                        ))
    }

    pub fn dressing_t_inv_landau<T: Num>(om: T, p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
        let sdim = om * om + p * p;
        let s = sdim / (m * m);
        gluon_vac::dressing_inv_landau(s, mq / m, f0)
            - sdim.inv()
                * PREFACTOR
                * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(om, p, m, beta)
                    + (nf() as R)
                        * gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                            om, p, mq, beta, mu,
                        ))
    }

    pub fn dressing_l_landau<T: Num>(om: T, p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
        dressing_l_inv_landau(om, p, m, mq, beta, mu, f0).inv()
    }

    pub fn dressing_t_landau<T: Num>(om: T, p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
        dressing_t_inv_landau(om, p, m, mq, beta, mu, f0).inv()
    }

    pub fn propagator_l_landau<T: Num>(om: T, p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
        (om * om + p * p).inv() * dressing_l_landau(om, p, m, mq, beta, mu, f0)
    }

    pub fn propagator_t_landau<T: Num>(om: T, p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
        (om * om + p * p).inv() * dressing_t_landau(om, p, m, mq, beta, mu, f0)
    }

    pub fn dressing_l_inv_zero_temp_landau<T: Num>(om: T, p: R, m: R, mq: R, mu: R, f0: R) -> C {
        let sdim = om * om + p * p;
        let s = sdim / (m * m);
        gluon_vac::dressing_inv_landau(s, mq / m, f0)
            - sdim.inv()
                * PREFACTOR
                * (nf() as R)
                * gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(
                    om, p, mq, mu,
                )
    }

    pub fn dressing_t_inv_zero_temp_landau<T: Num>(om: T, p: R, m: R, mq: R, mu: R, f0: R) -> C {
        let sdim = om * om + p * p;
        let s = sdim / (m * m);
        gluon_vac::dressing_inv_landau(s, mq / m, f0)
            - sdim.inv()
                * PREFACTOR
                * (nf() as R)
                * gluon_thermal_parts::polarization_quark_t_thermal_part_zero_temp_landau(
                    om, p, mq, mu,
                )
    }

    pub fn dressing_l_zero_temp_landau<T: Num>(om: T, p: R, m: R, mq: R, mu: R, f0: R) -> C {
        dressing_l_inv_zero_temp_landau(om, p, m, mq, mu, f0).inv()
    }

    pub fn dressing_t_zero_temp_landau<T: Num>(om: T, p: R, m: R, mq: R, mu: R, f0: R) -> C {
        dressing_t_inv_zero_temp_landau(om, p, m, mq, mu, f0).inv()
    }

    pub fn propagator_l_zero_temp_landau<T: Num>(om: T, p: R, m: R, mq: R, mu: R, f0: R) -> C {
        (om * om + p * p).inv() * dressing_l_zero_temp_landau(om, p, m, mq, mu, f0)
    }

    pub fn propagator_t_zero_temp_landau<T: Num>(om: T, p: R, m: R, mq: R, mu: R, f0: R) -> C {
        (om * om + p * p).inv() * dressing_t_zero_temp_landau(om, p, m, mq, mu, f0)
    }

    pub fn dressing_l_inv_landau_w_field_config<T1: Num, T2: Num>(
        om: T1,
        p: T2,
        beta: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        let m = config.gluon;
        let sdim = om * om + (p * p).into();
        let s = sdim / (m * m);
        let quark_contribution: C = config
            .quarks_internal
            .iter()
            .map(|qi| {
                gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                    om, p, qi.mq, beta, mu,
                ) * (qi.nf as R)
            })
            .sum();
        gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
            - sdim.inv()
                * PREFACTOR
                * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(om, p, m, beta)
                    + quark_contribution)
    }

    pub fn dressing_t_inv_landau_w_field_config<T: Num>(
        om: T,
        p: R,
        beta: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        let m = config.gluon;
        let sdim = om * om + p * p;
        let s = sdim / (m * m);
        let quark_contribution: C = config
            .quarks_internal
            .iter()
            .map(|qi| {
                gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                    om, p, qi.mq, beta, mu,
                ) * (qi.nf as R)
            })
            .sum();
        gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
            - sdim.inv()
                * PREFACTOR
                * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(om, p, m, beta)
                    + quark_contribution)
    }

    pub fn dressing_l_landau_w_field_config<T1: Num, T2: Num>(
        om: T1,
        p: T2,
        beta: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        dressing_l_inv_landau_w_field_config(om, p, beta, mu, f0, config).inv()
    }

    pub fn dressing_t_landau_w_field_config<T: Num>(
        om: T,
        p: R,
        beta: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        dressing_t_inv_landau_w_field_config(om, p, beta, mu, f0, config).inv()
    }

    pub fn propagator_l_landau_w_field_config<T1: Num, T2: Num>(
        om: T1,
        p: T2,
        beta: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        (om * om + (p * p).into()).inv()
            * dressing_l_landau_w_field_config(om, p, beta, mu, f0, config)
    }

    pub fn propagator_t_landau_w_field_config<T: Num>(
        om: T,
        p: R,
        beta: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        (om * om + p * p).inv() * dressing_t_landau_w_field_config(om, p, beta, mu, f0, config)
    }

    pub fn dressing_l_inv_zero_temp_landau_w_field_config<T1: Num, T2: Num>(
        om: T1,
        p: T2,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        let m = config.gluon;
        let sdim = om * om + (p * p).into();
        let s = sdim / (m * m);
        let quark_contribution: C = config
            .quarks_internal
            .iter()
            .map(|qi| {
                gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(
                    om, p, qi.mq, mu,
                ) * (qi.nf as R)
            })
            .sum();
        gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
            - sdim.inv() * PREFACTOR * quark_contribution
    }

    pub fn dressing_t_inv_zero_temp_landau_w_field_config<T: Num>(
        om: T,
        p: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        let m = config.gluon;
        let sdim = om * om + p * p;
        let s = sdim / (m * m);
        let quark_contribution: C = config
            .quarks_internal
            .iter()
            .map(|qi| {
                gluon_thermal_parts::polarization_quark_t_thermal_part_zero_temp_landau(
                    om, p, qi.mq, mu,
                ) * (qi.nf as R)
            })
            .sum();
        gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
            - sdim.inv() * PREFACTOR * quark_contribution
    }

    pub fn dressing_l_zero_temp_landau_w_field_config<T1: Num, T2: Num>(
        om: T1,
        p: T2,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        dressing_l_inv_zero_temp_landau_w_field_config(om, p, mu, f0, config).inv()
    }

    pub fn dressing_t_zero_temp_landau_w_field_config<T: Num>(
        om: T,
        p: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        dressing_t_inv_zero_temp_landau_w_field_config(om, p, mu, f0, config).inv()
    }

    pub fn propagator_l_zero_temp_landau_w_field_config<T1: Num, T2: Num>(
        om: T1,
        p: T2,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        (om * om + (p * p).into()).inv()
            * dressing_l_zero_temp_landau_w_field_config(om, p, mu, f0, config)
    }

    pub fn propagator_t_zero_temp_landau_w_field_config<T: Num>(
        om: T,
        p: R,
        mu: R,
        f0: R,
        config: &FieldConfig,
    ) -> C {
        (om * om + p * p).inv() * dressing_t_zero_temp_landau_w_field_config(om, p, mu, f0, config)
    }

    pub(crate) mod ffi {}
}

pub mod quark {
    use crate::common::thermal::inlines::fermi_momentum;
    use crate::consts::{get_default_integration_method, nc};
    use crate::utils::find_root;
    use crate::{types::Integral, R};
    use peroxide::numerical::integral::integrate;
    use std::f64::consts::PI;

    const PI2: R = PI * PI;

    pub fn charge_density_with_method(m: R, beta: R, mu: R, integral: Integral) -> R {
        integrate(
            |t| inlines::charge_density_i((1. - t) / t, m, beta, mu) / (t * t),
            (0., 1.),
            integral,
        )
    }

    pub fn charge_density(m: R, beta: R, mu: R) -> R {
        charge_density_with_method(m, beta, mu, get_default_integration_method())
    }

    pub fn charge_density_zero_temp(m: R, mu: R) -> R {
        if mu.abs() <= m.abs() {
            return 0.;
        }
        let s = if mu > m { 1. } else { -1. };
        let fm = fermi_momentum(m, mu);
        (nc() as R) * s * (fm * fm * fm) / (3. * PI2)
    }

    pub fn chemical_potential_with_method(m: R, beta: R, n: R, integral: Integral) -> R {
        if n == 0. {
            return 0.;
        }
        find_root(
            &|mu| n - charge_density_with_method(m, beta, mu, integral),
            (0.5 * n, n),
            integral.get_tol(),
            1000,
        )
        .unwrap()
        .0
    }

    pub fn chemical_potential(m: R, beta: R, n: R) -> R {
        chemical_potential_with_method(m, beta, n, get_default_integration_method())
    }

    pub fn chemical_potential_zero_temp(m: R, n: R) -> R {
        if n == 0. {
            return 0.;
        }
        let a = (3. * PI2 * n.abs() / (nc() as R)).powf(2. / 3.);
        let s = if n > 0. { 1. } else { -1. };
        s * m.mul_add(m, a).sqrt()
    }

    pub(crate) mod inlines {
        use crate::common::{inlines::energy, thermal::inlines::fermi_distribution};
        use crate::consts::nc;
        use crate::R;
        use std::f64::consts::PI;

        const PI2: R = PI * PI;

        pub fn charge_density_i(q: R, m: R, beta: R, mu: R) -> R {
            let en = energy(q, m);
            (nc() as R)
                * q
                * q
                * (fermi_distribution(en, beta, mu) - fermi_distribution(en, beta, -mu))
                / (PI2)
        }
    }
}

pub mod zero_matsubara {
    pub mod ghost {
        pub use crate::ym::thermal::zero_matsubara::ghost::*;
    }

    pub mod gluon {
        use crate::consts::nf;
        use crate::low_level::oneloop::thermal::gluon::zero_matsubara as gluon_thermal_parts;
        use crate::qcd::gluon as gluon_vac;
        use crate::qcd::FieldConfig;
        use crate::R;
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        pub fn dressing_l_inv_landau(p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> R {
            let sdim = p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - PREFACTOR
                    * (nf() as R).mul_add(
                        gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                            p, mq, beta, mu,
                        ),
                        gluon_thermal_parts::polarization_glue_l_thermal_part_landau(p, m, beta),
                    )
                    / sdim
        }

        pub fn dressing_t_inv_landau(p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> R {
            let sdim = p * p;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - PREFACTOR
                    * (nf() as R).mul_add(
                        gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                            p, mq, beta, mu,
                        ),
                        gluon_thermal_parts::polarization_glue_t_thermal_part_landau(p, m, beta),
                    )
                    / sdim
        }

        pub fn dressing_l_landau(p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> R {
            1. / dressing_l_inv_landau(p, m, mq, beta, mu, f0)
        }

        pub fn dressing_t_landau(p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> R {
            1. / dressing_t_inv_landau(p, m, mq, beta, mu, f0)
        }

        pub fn propagator_l_landau(p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> R {
            dressing_l_landau(p, m, mq, beta, mu, f0) / (p * p)
        }

        pub fn propagator_t_landau(p: R, m: R, mq: R, beta: R, mu: R, f0: R) -> R {
            dressing_t_landau(p, m, mq, beta, mu, f0) / (p * p)
        }

        pub fn dressing_l_inv_zero_temp_landau(p: R, m: R, mq: R, mu: R, f0: R) -> R {
            let s = p * p / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - PREFACTOR
                    * (nf() as R)
                    * gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(
                        p, mq, mu,
                    )
                    / (p * p)
        }

        pub fn dressing_t_inv_zero_temp_landau(p: R, m: R, mq: R, mu: R, f0: R) -> R {
            let s = p * p / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - PREFACTOR
                    * (nf() as R)
                    * gluon_thermal_parts::polarization_quark_t_thermal_part_zero_temp_landau(
                        p, mq, mu,
                    )
                    / (p * p)
        }

        pub fn dressing_l_zero_temp_landau(p: R, m: R, mq: R, mu: R, f0: R) -> R {
            1. / dressing_l_inv_zero_temp_landau(p, m, mq, mu, f0)
        }

        pub fn dressing_t_zero_temp_landau(p: R, m: R, mq: R, mu: R, f0: R) -> R {
            1. / dressing_t_inv_zero_temp_landau(p, m, mq, mu, f0)
        }

        pub fn propagator_l_zero_temp_landau(p: R, m: R, mq: R, mu: R, f0: R) -> R {
            dressing_l_zero_temp_landau(p, m, mq, mu, f0) / (p * p)
        }

        pub fn propagator_t_zero_temp_landau(p: R, m: R, mq: R, mu: R, f0: R) -> R {
            dressing_t_zero_temp_landau(p, m, mq, mu, f0) / (p * p)
        }

        pub fn dressing_l_inv_landau_w_field_config(
            p: R,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            let m = config.gluon;
            let sdim = p * p;
            let s = sdim / (m * m);
            let quark_contribution: R = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                        p, qi.mq, beta, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(p, m, beta)
                        + quark_contribution)
                    / sdim
        }

        pub fn dressing_t_inv_landau_w_field_config(
            p: R,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            let m = config.gluon;
            let sdim = p * p;
            let s = sdim / (m * m);
            let quark_contribution: R = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                        p, qi.mq, beta, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(p, m, beta)
                        + quark_contribution)
                    / sdim
        }

        pub fn dressing_l_landau_w_field_config(
            p: R,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            1. / dressing_l_inv_landau_w_field_config(p, beta, mu, f0, config)
        }

        pub fn dressing_t_landau_w_field_config(
            p: R,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            1. / dressing_t_inv_landau_w_field_config(p, beta, mu, f0, config)
        }

        pub fn propagator_l_landau_w_field_config(
            p: R,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            dressing_l_landau_w_field_config(p, beta, mu, f0, config) / (p * p)
        }

        pub fn propagator_t_landau_w_field_config(
            p: R,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            dressing_t_landau_w_field_config(p, beta, mu, f0, config) / (p * p)
        }

        pub fn dressing_l_inv_zero_temp_landau_w_field_config(
            p: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            let m = config.gluon;
            let s = p * p / (m * m);
            let quark_contribution: R = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(
                        p, qi.mq, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - PREFACTOR * quark_contribution / (p * p)
        }

        pub fn dressing_t_inv_zero_temp_landau_w_field_config(
            p: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            let m = config.gluon;
            let s = p * p / (m * m);
            let quark_contribution: R = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_t_thermal_part_zero_temp_landau(
                        p, qi.mq, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - PREFACTOR * quark_contribution / (p * p)
        }

        pub fn dressing_l_zero_temp_landau_w_field_config(
            p: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            1. / dressing_l_inv_zero_temp_landau_w_field_config(p, mu, f0, config)
        }

        pub fn dressing_t_zero_temp_landau_w_field_config(
            p: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            1. / dressing_t_inv_zero_temp_landau_w_field_config(p, mu, f0, config)
        }

        pub fn propagator_l_zero_temp_landau_w_field_config(
            p: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            dressing_l_zero_temp_landau_w_field_config(p, mu, f0, config) / (p * p)
        }

        pub fn propagator_t_zero_temp_landau_w_field_config(
            p: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> R {
            dressing_t_zero_temp_landau_w_field_config(p, mu, f0, config) / (p * p)
        }

        pub(crate) mod ffi {}
    }
}

pub mod zero_momentum {
    pub mod ghost {
        pub use crate::ym::thermal::zero_momentum::ghost::*;
    }

    pub mod gluon {
        use crate::consts::nf;
        use crate::low_level::oneloop::thermal::gluon::zero_momentum as gluon_thermal_parts;
        use crate::qcd::gluon as gluon_vac;
        use crate::qcd::FieldConfig;
        use crate::{Num, C, R};
        use std::f64::consts::PI;

        const PREFACTOR: R = (16. * PI * PI) / 3.;

        pub fn dressing_l_inv_landau<T: Num>(om: T, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
            let sdim = om * om;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - sdim.inv()
                    * PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(om, m, beta)
                        + (nf() as R)
                            * gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                                om, mq, beta, mu,
                            ))
        }

        pub fn dressing_t_inv_landau<T: Num>(om: T, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
            let sdim = om * om;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - sdim.inv()
                    * PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(om, m, beta)
                        + (nf() as R)
                            * gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                                om, mq, beta, mu,
                            ))
        }

        pub fn dressing_l_landau<T: Num>(om: T, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
            dressing_l_inv_landau(om, m, mq, beta, mu, f0).inv()
        }

        pub fn dressing_t_landau<T: Num>(om: T, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
            dressing_t_inv_landau(om, m, mq, beta, mu, f0).inv()
        }

        pub fn propagator_l_landau<T: Num>(om: T, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
            (om * om).inv() * dressing_l_landau(om, m, mq, beta, mu, f0)
        }

        pub fn propagator_t_landau<T: Num>(om: T, m: R, mq: R, beta: R, mu: R, f0: R) -> C {
            (om * om).inv() * dressing_t_landau(om, m, mq, beta, mu, f0)
        }

        pub fn dressing_l_inv_zero_temp_landau<T: Num>(om: T, m: R, mq: R, mu: R, f0: R) -> C {
            let sdim = om * om;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - sdim.inv()
                    * PREFACTOR
                    * (nf() as R)
                    * gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(
                        om, mq, mu,
                    )
        }

        pub fn dressing_t_inv_zero_temp_landau<T: Num>(om: T, m: R, mq: R, mu: R, f0: R) -> C {
            let sdim = om * om;
            let s = sdim / (m * m);
            gluon_vac::dressing_inv_landau(s, mq / m, f0)
                - sdim.inv()
                    * PREFACTOR
                    * (nf() as R)
                    * gluon_thermal_parts::polarization_quark_t_thermal_part_zero_temp_landau(
                        om, mq, mu,
                    )
        }

        pub fn dressing_l_zero_temp_landau<T: Num>(om: T, m: R, mq: R, mu: R, f0: R) -> C {
            dressing_l_inv_zero_temp_landau(om, m, mq, mu, f0).inv()
        }

        pub fn dressing_t_zero_temp_landau<T: Num>(om: T, m: R, mq: R, mu: R, f0: R) -> C {
            dressing_t_inv_zero_temp_landau(om, m, mq, mu, f0).inv()
        }

        pub fn propagator_l_zero_temp_landau<T: Num>(om: T, m: R, mq: R, mu: R, f0: R) -> C {
            (om * om).inv() * dressing_l_zero_temp_landau(om, m, mq, mu, f0)
        }

        pub fn propagator_t_zero_temp_landau<T: Num>(om: T, m: R, mq: R, mu: R, f0: R) -> C {
            (om * om).inv() * dressing_t_zero_temp_landau(om, m, mq, mu, f0)
        }

        pub fn dressing_l_inv_landau_w_field_config<T: Num>(
            om: T,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            let m = config.gluon;
            let sdim = om * om;
            let s = sdim / (m * m);
            let quark_contribution: C = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_l_thermal_part_landau(
                        om, qi.mq, beta, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - sdim.inv()
                    * PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_l_thermal_part_landau(om, m, beta)
                        + quark_contribution)
        }

        pub fn dressing_t_inv_landau_w_field_config<T: Num>(
            om: T,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            let m = config.gluon;
            let sdim = om * om;
            let s = sdim / (m * m);
            let quark_contribution: C = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_t_thermal_part_landau(
                        om, qi.mq, beta, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - sdim.inv()
                    * PREFACTOR
                    * (gluon_thermal_parts::polarization_glue_t_thermal_part_landau(om, m, beta)
                        + quark_contribution)
        }

        pub fn dressing_l_landau_w_field_config<T: Num>(
            om: T,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            dressing_l_inv_landau_w_field_config(om, beta, mu, f0, config).inv()
        }

        pub fn dressing_t_landau_w_field_config<T: Num>(
            om: T,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            dressing_t_inv_landau_w_field_config(om, beta, mu, f0, config).inv()
        }

        pub fn propagator_l_landau_w_field_config<T: Num>(
            om: T,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            (om * om).inv() * dressing_l_landau_w_field_config(om, beta, mu, f0, config)
        }

        pub fn propagator_t_landau_w_field_config<T: Num>(
            om: T,
            beta: R,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            (om * om).inv() * dressing_t_landau_w_field_config(om, beta, mu, f0, config)
        }

        pub fn dressing_l_inv_zero_temp_landau_w_field_config<T: Num>(
            om: T,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            let m = config.gluon;
            let sdim = om * om;
            let s = sdim / (m * m);
            let quark_contribution: C = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(
                        om, qi.mq, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - sdim.inv() * PREFACTOR * quark_contribution
        }

        pub fn dressing_t_inv_zero_temp_landau_w_field_config<T: Num>(
            om: T,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            let m = config.gluon;
            let sdim = om * om;
            let s = sdim / (m * m);
            let quark_contribution: C = config
                .quarks_internal
                .iter()
                .map(|qi| {
                    gluon_thermal_parts::polarization_quark_t_thermal_part_zero_temp_landau(
                        om, qi.mq, mu,
                    ) * (qi.nf as R)
                })
                .sum();
            gluon_vac::dressing_inv_landau_w_field_config(s, f0, config)
                - sdim.inv() * PREFACTOR * quark_contribution
        }

        pub fn dressing_l_zero_temp_landau_w_field_config<T: Num>(
            om: T,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            dressing_l_inv_zero_temp_landau_w_field_config(om, mu, f0, config).inv()
        }

        pub fn dressing_t_zero_temp_landau_w_field_config<T: Num>(
            om: T,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            dressing_t_inv_zero_temp_landau_w_field_config(om, mu, f0, config).inv()
        }

        pub fn propagator_l_zero_temp_landau_w_field_config<T: Num>(
            om: T,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            (om * om).inv() * dressing_l_zero_temp_landau_w_field_config(om, mu, f0, config)
        }

        pub fn propagator_t_zero_temp_landau_w_field_config<T: Num>(
            om: T,
            mu: R,
            f0: R,
            config: &FieldConfig,
        ) -> C {
            (om * om).inv() * dressing_t_zero_temp_landau_w_field_config(om, mu, f0, config)
        }

        pub(crate) mod ffi {}
    }
}

#[cfg(test)]
mod tests {
    use crate::{Num, R};

    const TOLERANCE: R = 10E-12;

    fn assert_equal_with_tol<T: Num>(lhs: T, rhs: T, tol: R) {
        if rhs != T::zero() {
            assert!(
                (lhs / rhs - 1.).abs() < tol,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {tol:e} rhs"
            );
        } else {
            assert!(
                (lhs - rhs).abs() < tol,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {tol:e}"
            );
        }
    }

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
    fn test_chemical_potential() {
        use crate::consts::get_default_quark_mass;
        use crate::qcd::thermal::quark::{
            charge_density_with_method, charge_density_zero_temp, chemical_potential_with_method,
            chemical_potential_zero_temp,
        };
        use crate::types::Integral;

        // Decrease tolerance to speed up tests. Do not decrease it too much though.
        let tol = 1E-8;
        let iter = 20;
        let m = get_default_quark_mass();

        for mu in [0.8, 2.5, 3.8] {
            for beta in [0.1, 0.2, 0.5] {
                assert_equal_with_tol(
                    chemical_potential_with_method(
                        m,
                        beta,
                        charge_density_with_method(m, beta, mu, Integral::G7K15(tol, iter)),
                        Integral::G7K15(tol, iter),
                    ),
                    mu,
                    tol,
                );
            }

            assert_equal(
                chemical_potential_zero_temp(m, charge_density_zero_temp(m, mu)),
                mu,
            );
        }
    }
}
