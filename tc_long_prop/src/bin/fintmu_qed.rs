use gnuplot::{AutoOption, AxesCommon, Figure, LabelOption};
use lazy_static::lazy_static;
use log::info;
#[cfg(feature = "ftmq_alt_riemann")]
use qcd_sme::common::thermal::{
    fermi_distribution_double as fermi_distribution,
    fermi_distribution_double_zero_temp as fermi_distribution_zero_temp,
};
use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::low_level::oneloop::gluon::f_q;
use qcd_sme::low_level::oneloop::thermal::gluon as gluon_thermal_parts;
use qcd_sme::{qcd::FieldConfig, Num, C, R};
use rayon::prelude::*;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;
use tc_long_prop::BASEDIR;

const PREFACTOR: R = (16. * PI * PI) / 3.;

// This is not really used, it cancels out from expressions
const MG: R = 1.;
const F0: R = 2.89;

const P0: R = 0.01;
const PREN: R = 4.;
const OMEPS: R = 1E-3;
const OMEPSC: C = C::new(OMEPS, 0.);

const PIXX: u32 = 1600;
const PIXY: u32 = 1600;

const SHOW: bool = false;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_qed");
}

fn max_idx(mat: &[R]) -> usize {
    let (mut mx, mut imx) = (0., 0);
    for (i, &m) in mat.iter().enumerate() {
        if m > mx {
            mx = m;
            imx = i;
        }
    }
    imx
}

// Adapted from qcd_sme_rs::qcd::gluon::inlines::dressing_inv_landau_w_field_config
fn qed_dressing_vac_inv(s: C, f0: R, config: &FieldConfig) -> C {
    let nc = config.nc;
    let mg = config.gluon;
    let mg2 = mg * mg;
    let mut res = C::new(f0, 0.);
    for (nf, mq) in &config.quarks {
        let s_q = s / (mq * mq / mg2);
        res += f_q(s_q) * (*nf as R / nc as R);
    }
    res
}

#[cfg(feature = "ftmq_alt_riemann")]
fn qed_riemann_corr_zero_temp(om: C, mu: R) -> C {
    C::I * om * om / (12. * PI) * (fermi_distribution_zero_temp(C::I * om / 2., mu) - 1.) / 3.
}

#[cfg(feature = "ftmq_alt_riemann")]
fn qed_riemann_corr(om: C, beta: R, mu: R) -> C {
    C::I * om * om / (12. * PI) * (fermi_distribution(C::I * om / 2., beta, mu) - 1.) / 3.
}

// Adapted from qcd_sme_rs::qcd::thermal::gluon::dressing_t_inv_zero_temp_landau_w_field_config
fn qed_dressing_t_zero_temp(om: C, p: R, mu: R, f0: R, config: &FieldConfig) -> C {
    let m = config.gluon;
    let sdim = om * om + p * p;
    let s = sdim / (m * m);
    #[allow(unused_mut)]
    let mut quark_contribution: C = config
        .quarks
        .iter()
        .map(|&(nf, mq)| {
            gluon_thermal_parts::polarization_quark_l_thermal_part_zero_temp_landau(om, p, mq, mu)
                * (nf as R)
        })
        .sum();
    #[cfg(feature = "ftmq_alt_riemann")]
    if om.re < 0. {
        quark_contribution += qed_riemann_corr_zero_temp(om, mu);
    }
    1. / (qed_dressing_vac_inv(s, f0, config) - sdim.inv() * PREFACTOR * quark_contribution)
}

fn qed_prop_t_zero_temp(om: C, p: R, mu: R, f0: R, config: &FieldConfig) -> C {
    let sdim = om * om + p * p;
    qed_dressing_t_zero_temp(om, p, mu, f0, config) / sdim
}

// Adapted from qcd_sme_rs::qcd::thermal::gluon::dressing_t_inv_landau_w_field_config
fn qed_dressing_t(om: C, p: R, beta: R, mu: R, f0: R, config: &FieldConfig) -> C {
    let m = config.gluon;
    let sdim = om * om + p * p;
    let s = sdim / (m * m);
    #[allow(unused_mut)]
    let mut quark_contribution: C = config
        .quarks
        .iter()
        .map(|&(nf, mq)| {
            gluon_thermal_parts::polarization_quark_t_thermal_part_landau(om, p, mq, beta, mu)
                * (nf as R)
        })
        .sum();
    #[cfg(feature = "ftmq_alt_riemann")]
    if om.re < 0. {
        quark_contribution += qed_riemann_corr(om, beta, mu);
    }
    1. / (qed_dressing_vac_inv(s, f0, config) - sdim.inv() * PREFACTOR * quark_contribution)
}

fn qed_prop_t(om: C, p: R, beta: R, mu: R, f0: R, config: &FieldConfig) -> C {
    let sdim = om * om + p * p;
    qed_dressing_t(om, p, beta, mu, f0, config) / sdim
}

// Adapted from fintmu_spectral::plot_qcd_t
fn plot_qed_spectral_t(oms: &[R], t: R, config: &FieldConfig, f0: R, dir: &str) -> R {
    let spectr: Vec<R> = if t == 0. {
        let z = 1. / qed_dressing_t_zero_temp(OMEPSC, PREN, 0., f0, config).re;
        oms.par_iter()
            .map(| &om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*qed_prop_t_zero_temp(om, P0, 0., f0, config).im/PI
                };
                info!(
                    "Computed QED mu=0 spectral function at T = 0.000 GeV for om = {om:.4} GeV: {res}"
                );
                res
            })
            .collect()
    } else {
        let z = 1. / qed_dressing_t(OMEPSC, PREN, 1. / t, 0., f0, config).re;
        oms.par_iter()
            .map(|&om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*qed_prop_t(om, P0, 1. / t, 0., f0, config).im/PI
                };
                info!(
                    "Computed QED mu=0 spectral function at T = {t:.3} GeV for om = {om:.4} GeV: {res}"
                );
                res
            })
            .collect()
    };

    let mut figure = Figure::new();
    figure
        .axes2d()
        .lines(oms, &spectr, &[])
        .set_x_range(
            AutoOption::Fix(*oms.first().unwrap()),
            AutoOption::Fix(*oms.last().unwrap()),
        )
        .set_y_range(AutoOption::Fix(-0.5), AutoOption::Fix(50.));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qed_spectr_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    oms[max_idx(&spectr)]
}

fn plot_qed_dressing_t(momenta: &[R], t: R, config: &FieldConfig, f0: R, dir: &str) {
    let prop: Vec<R> = if t == 0. {
        let z = 1. / qed_dressing_t_zero_temp(OMEPSC, PREN, 0., f0, config).re;
        momenta
            .par_iter()
            .map(|&p| {
                let res = if p == 0. {
                    R::NAN
                } else {
                    z * qed_dressing_t_zero_temp(OMEPSC, p, 0., f0, config).re
                };
                info!(
                    "Computed QED mu=0 dressing function at T = 0.000 GeV for p = {p:.4} GeV: {res}"
                );
                res
            })
            .collect()
    } else {
        let z = 1. / qed_dressing_t(OMEPSC, PREN, 1. / t, 0., f0, config).re;
        momenta
            .par_iter()
            .map(|&p| {
                let res = if p == 0. {
                    R::NAN
                } else {
                    z * qed_dressing_t(OMEPSC, p, 1. / t, 0., f0, config).re
                };
                info!(
                    "Computed QED mu=0 dressing function at T = {t:.3} GeV for p = {p:.4} GeV: {res}"
                );
                res
            })
            .collect()
    };

    let mut figure = Figure::new();
    figure
        .axes2d()
        .lines(momenta, &prop, &[])
        .set_x_range(
            AutoOption::Fix(*momenta.first().unwrap()),
            AutoOption::Fix(*momenta.last().unwrap()),
        )
        .set_y_range(AutoOption::Fix(0.), AutoOption::Fix(2.));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qed_dress_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }
}

#[cfg(feature = "ftmq_3d_poles")]
#[allow(clippy::too_many_arguments)]
fn plot_qed_dressing_cplane_t(
    oms: &[C],
    t: R,
    f0: R,
    config: &FieldConfig,
    ommax: R,
    nom: usize,
    omimscale: R,
    dir: &str,
) -> C {
    let nrc = nom + 1;
    let nact = nrc * nrc;

    let mat: Vec<R> = if t == 0. {
        oms.par_iter().enumerate()
            .map(|(i, &om)| {
                let res = if om.re == 0. { R::NAN } else {
                    qed_dressing_t_zero_temp(om, P0, 0., f0, config).abs()
                };
                info!(
                    "Computed QED dressing function at T = 0.000 GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1
                );
                res
            })
            .collect()
    } else {
        oms.par_iter().enumerate()
            .map(|(i,&om)| {
                let res = if om.re == 0. { R::NAN } else {
                    qed_dressing_t(om, P0,1./t, 0., f0, config).abs()
                };
                info!(
                    "Computed QED dressing function at T = {t:.3} GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1
                );
                res
            })
            .collect()
    };

    let omx = oms[max_idx(&mat)];

    #[cfg(feature = "ftmq_alt_riemann")]
    let (zmin, zmax) = (0., 50.);
    #[cfg(not(feature = "ftmq_alt_riemann"))]
    let (zmin, zmax) = (0., 20.);

    let mut figure = Figure::new();
    figure
        .axes3d()
        .surface(
            mat,
            nrc,
            nrc,
            Some((-ommax / omimscale, -ommax, ommax / omimscale, ommax)),
            &[],
        )
        .set_x_range(
            AutoOption::Fix(-ommax / omimscale),
            AutoOption::Fix(ommax / omimscale),
        )
        .set_y_range(AutoOption::Fix(-ommax), AutoOption::Fix(ommax))
        .set_z_range(AutoOption::Fix(zmin), AutoOption::Fix(zmax))
        .set_x_label("Im(\u{03c9})", &[])
        .set_y_label("Re(\u{03c9})", &[])
        .set_title(
            &format!(
                "T = {t:.3} GeV, \u{03c9} = ({z:.4}) GeV",
                z = C::new(omx.im, omx.re)
            ),
            &[
                LabelOption::Font("Arial", 30.),
                LabelOption::TextOffset(0., -10.),
            ],
        );

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qed_dress_cplane_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    omx
}

fn plot_qed_spectral_mu(oms: &[R], mu: R, config: &FieldConfig, f0: R, dir: &str) -> R {
    let z = 1. / qed_dressing_t_zero_temp(OMEPSC, PREN, mu, f0, config).re;
    let spectr: Vec<R> = oms.par_iter()
            .map(| &om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*qed_prop_t_zero_temp(om, P0, mu, f0, config).im/PI
                };
                info!(
                    "Computed QED T=0 spectral function at mu = {mu:.3} GeV for om = {om:.4} GeV: {res}"
                );
                res
            })
            .collect();

    let mut figure = Figure::new();
    figure
        .axes2d()
        .lines(oms, &spectr, &[])
        .set_x_range(
            AutoOption::Fix(*oms.first().unwrap()),
            AutoOption::Fix(*oms.last().unwrap()),
        )
        .set_y_range(AutoOption::Fix(-0.5), AutoOption::Fix(50.));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qed_spectr_mu{mu:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    oms[max_idx(&spectr)]
}

fn plot_qed_dressing_mu(momenta: &[R], mu: R, config: &FieldConfig, f0: R, dir: &str) {
    let z = 1. / qed_dressing_t_zero_temp(OMEPSC, PREN, mu, f0, config).re;
    let prop: Vec<R> = momenta
        .par_iter()
        .map(|&p| {
            let res = if p == 0. {
                R::NAN
            } else {
                z * qed_dressing_t_zero_temp(OMEPSC, p, mu, f0, config).re
            };
            info!(
                "Computed QED T=0 dressing function at mu = {mu:.3} GeV for p = {p:.4} GeV: {res}"
            );
            res
        })
        .collect();

    let mut figure = Figure::new();
    figure
        .axes2d()
        .lines(momenta, &prop, &[])
        .set_x_range(
            AutoOption::Fix(*momenta.first().unwrap()),
            AutoOption::Fix(*momenta.last().unwrap()),
        )
        .set_y_range(AutoOption::Fix(0.), AutoOption::Fix(2.));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qed_dress_mu{mu:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }
}

#[cfg(feature = "ftmq_3d_poles")]
#[allow(clippy::too_many_arguments)]
fn plot_qed_dressing_cplane_mu(
    oms: &[C],
    mu: R,
    f0: R,
    config: &FieldConfig,
    ommax: R,
    nom: usize,
    omimscale: R,
    dir: &str,
) -> C {
    let nrc = nom + 1;
    let nact = nrc * nrc;

    let mat: Vec<R> =
        oms.par_iter().enumerate()
            .map(|(i, &om)| {
                let res = if om.re == 0. { R::NAN } else {
                    qed_dressing_t_zero_temp(om, P0, mu, f0, config).abs()
                };
                info!(
                    "Computed QED dressing function at mu = {mu:.3} GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1
                );
                res
            })
            .collect()
    ;

    let omx = oms[max_idx(&mat)];

    #[cfg(feature = "ftmq_alt_riemann")]
    let (zmin, zmax) = (0., 50.);
    #[cfg(not(feature = "ftmq_alt_riemann"))]
    let (zmin, zmax) = (0., 20.);

    let mut figure = Figure::new();
    figure
        .axes3d()
        .surface(
            mat,
            nrc,
            nrc,
            Some((-ommax / omimscale, -ommax, ommax / omimscale, ommax)),
            &[],
        )
        .set_x_range(
            AutoOption::Fix(-ommax / omimscale),
            AutoOption::Fix(ommax / omimscale),
        )
        .set_y_range(AutoOption::Fix(-ommax), AutoOption::Fix(ommax))
        .set_z_range(AutoOption::Fix(zmin), AutoOption::Fix(zmax))
        .set_x_label("Im(\u{03c9})", &[])
        .set_y_label("Re(\u{03c9})", &[])
        .set_title(
            &format!(
                "mu = {mu:.3} GeV, \u{03c9} = ({z:.4}) GeV",
                z = C::new(omx.im, omx.re)
            ),
            &[
                LabelOption::Font("Arial", 30.),
                LabelOption::TextOffset(0., -10.),
            ],
        );

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qed_dress_cplane_mu{mu:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    omx
}

fn qed_alpha(f0: R, t: R, config: &FieldConfig) -> R {
    if t == 0. {
        PREFACTOR / 24. / PI * PREN * PREN * qed_prop_t_zero_temp(OMEPSC, PREN, 0., f0, config).re
    } else {
        PREFACTOR / 24. / PI * PREN * PREN * qed_prop_t(OMEPSC, PREN, 1. / t, 0., f0, config).re
    }
}

fn qed_pole_literature(f0: R, t: R, config: &FieldConfig) -> R {
    (4. * PI / 9. * qed_alpha(f0, t, config)).sqrt() * t
}

fn qed_width_literature(f0: R, t: R, config: &FieldConfig) -> R {
    qed_alpha(f0, t, config) / 6. * qed_pole_literature(f0, t, config)
}

fn main() {
    /* 0. Setup */
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    set_default_max_iter_integral(20);

    if !SHOW {
        fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
        fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();
    };

    let ommin = 0.1;
    let ommax = 2.0;
    let nom = 200;
    let dom = (ommax - ommin) / (nom as R);

    let oms: Vec<R> = (0..=nom).map(|i| ommin + (i as R) * dom).collect();

    let pmin = 0.05;
    let pmax = 2.0;
    let np = 200;
    let dp = (pmax - pmin) / (np as R);

    let momenta: Vec<R> = (0..=np).map(|i| pmin + (i as R) * dp).collect();

    let omcmax = 1.0;
    let omcimscale = 50.;
    let nomc = 200;
    let domc = 2. * omcmax / (nomc as R);
    let nomc_actual = nomc * nomc + 2 * nomc + 1;

    let mut omcs = Vec::with_capacity(nomc_actual);

    for i in 0..=nomc {
        for j in 0..=nomc {
            omcs.push(C::new(
                ((j as R) * domc - omcmax) / omcimscale,
                (i as R) * domc - omcmax,
            ));
        }
    }

    // One single light quark. nc = 3 is needed for consistency elsewhere
    let config = FieldConfig::new(3, MG, vec![(1, 0.001)]);

    /* I. QED */

    // IA. mu = 0, T-dependent, fixed parameters
    if !SHOW {
        fs::create_dir_all(THIS_BASEDIR.join("qed_fixed_zero_density").as_path()).unwrap();
    }

    let (mut qed_alphas, mut qed_poles_lit, mut qed_widths_lit, mut qed_poles_comp) =
        (vec![], vec![], vec![], vec![]);

    [
        0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
        0.85, 0.9, 0.95, 1.,
    ]
    .iter()
    .for_each(|&t| {
        qed_poles_comp.push(plot_qed_spectral_t(
            &oms,
            t,
            &config,
            F0,
            "qed_fixed_zero_density",
        ));
        plot_qed_dressing_t(&momenta, t, &config, F0, "qed_fixed_zero_density");
        plot_qed_dressing_cplane_t(
            &omcs,
            t,
            F0,
            &config,
            omcmax,
            nomc,
            omcimscale,
            "qed_fixed_zero_density",
        );
        qed_alphas.push(qed_alpha(F0, t, &config));
        qed_poles_lit.push(qed_pole_literature(F0, t, &config));
        qed_widths_lit.push(qed_width_literature(F0, t, &config));
    });

    info!("Computed for alphas: {qed_alphas:.5?}");
    info!("Corresponding poles (computed): {qed_poles_comp:.5?}");
    info!("Corresponding poles (literature): {qed_poles_lit:.5?}");
    info!("Corresponding widths (literature): {qed_widths_lit:.5?}");

    // IB. T = 0, mu-dependent, fixed parameters
    if !SHOW {
        fs::create_dir_all(THIS_BASEDIR.join("qed_fixed_zero_temp").as_path()).unwrap();
    }

    let mut qed_poles_comp = vec![];

    [0., 0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
        .iter()
        .for_each(|&mu| {
            qed_poles_comp.push(plot_qed_spectral_mu(
                &oms,
                mu,
                &config,
                F0,
                "qed_fixed_zero_temp",
            ));
            plot_qed_dressing_mu(&momenta, mu, &config, F0, "qed_fixed_zero_temp");
            plot_qed_dressing_cplane_mu(
                &omcs,
                mu,
                F0,
                &config,
                omcmax,
                nomc,
                omcimscale,
                "qed_fixed_zero_temp",
            );
        });
}
