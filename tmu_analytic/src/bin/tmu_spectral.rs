use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::qcd::thermal::gluon::{
    propagator_t_landau_w_field_config as qcd_propagator,
    propagator_t_zero_temp_landau_w_field_config as qcd_propagator_zero_temp,
};
use qcd_sme::qcd::FieldConfig;
use qcd_sme::types::{C, NCTYPE, R};
use qcd_sme::ym::gluon::propagator_landau as ym_propagator_zero_temp;
use qcd_sme::ym::thermal::gluon::propagator_t_landau as ym_propagator;
use tmu_analytic::BASEDIR;

use gnuplot::{AutoOption, AxesCommon, Figure};
use lazy_static::lazy_static;
use log::info;
use rayon::prelude::*;
use std::f64::consts::PI;
use std::io::{self, Write};
use std::{fs, path::PathBuf};

const NC: NCTYPE = 3;
const MG: R = 0.656;
const MQ: R = 0.4;
const P0: R = 0.01;
const F0: R = -0.876;
const PREN: R = 4.;
const OMEPS: R = 1E-3;

const PIXX: u32 = 1600;
const PIXY: u32 = 1600;

const NO_QCD: bool = false;
const NO_YM: bool = false;
const SHOW: bool = false;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_spectral");
}

fn plot_ym_t(oms: &[R], t: R, m: R, f0: R, dir: &str) {
    let spectr: Vec<R> = if t == 0. {
        let z = 1.
            / ym_propagator_zero_temp((OMEPS * OMEPS + PREN * PREN) / (m * m), f0)
            / (PREN * PREN);
        oms.par_iter()
            .map(| &om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*ym_propagator_zero_temp((om * om + P0 * P0) / (m * m), f0).im/PI
                };
                info!(
                    "Computed pure Yang-Mills spectral function at T = 0.000 GeV for om = {om:.4} GeV: {res}"
                );
                res
            })
            .collect()
    } else {
        let z = 1. / ym_propagator(OMEPS, PREN, m, 1. / t, f0).re / (PREN * PREN);
        oms.par_iter()
            .map(|&om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*ym_propagator(om, P0, m, 1. / t, f0).im/PI
                };
                info!(
                    "Computed pure Yang-Mills spectral function at T = {t:.3} GeV for om = {om:.4} GeV: {res}"
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
        .set_y_range(AutoOption::Fix(-1.1), AutoOption::Fix(0.5));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("ym_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();

        let mut fout = io::BufWriter::new(
            fs::File::create(THIS_BASEDIR.join(dir).join(format!("ym_t{t:.3}.out"))).unwrap(),
        );
        oms.iter().zip(spectr.iter()).for_each(|(om, s)| {
            writeln!(fout, "{om:.4}\t{s:.4}").unwrap();
        });
    }
}

fn plot_qcd_t(oms: &[R], t: R, config: &FieldConfig, f0: R, dir: &str) {
    let spectr: Vec<R> = if t == 0. {
        let z = 1. / qcd_propagator_zero_temp(OMEPS, PREN, 0., f0, config).re / (PREN * PREN);
        oms.par_iter()
            .map(| &om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*qcd_propagator_zero_temp(om, P0, 0., f0, config).im/PI
                };
                info!(
                    "Computed full QCD mu=0 spectral function at T = 0.000 GeV for om = {om:.4} GeV: {res}"
                );
                res
            })
            .collect()
    } else {
        let z = 1. / qcd_propagator(OMEPS, PREN, 1. / t, 0., f0, config).re / (PREN * PREN);
        oms.par_iter()
            .map(|&om| {
                let res = if om == 0. { R::NAN } else {
                    let om = C::new(OMEPS, -om);
                    z*qcd_propagator(om, P0, 1. / t, 0., f0, config).im/PI
                };
                info!(
                    "Computed full QCD mu=0 spectral function at T = {t:.3} GeV for om = {om:.4} GeV: {res}"
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
        .set_y_range(AutoOption::Fix(-1.1), AutoOption::Fix(2.5));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qcd_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();

        let mut fout = io::BufWriter::new(
            fs::File::create(THIS_BASEDIR.join(dir).join(format!("qcd_t{t:.3}.out"))).unwrap(),
        );
        oms.iter().zip(spectr.iter()).for_each(|(om, s)| {
            writeln!(fout, "{om:.4}\t{s:.4}").unwrap();
        });
    }
}

fn plot_qcd_mu(oms: &[R], mu: R, config: &FieldConfig, f0: R, dir: &str) {
    let z = 1. / qcd_propagator_zero_temp(OMEPS, PREN, mu, f0, config).re / (PREN * PREN);
    let spectr: Vec<R> = oms
        .par_iter()
        .map(|&om| {
            let res = if om == 0. {
                R::NAN
            } else {
                let om = C::new(OMEPS, -om);
                z * qcd_propagator_zero_temp(om, P0, mu, f0, config).im / PI
            };
            info!(
                "Computed full QCD T=0 spectral function at mu = {mu:.3} GeV for om = {om:.4} GeV: {res}"
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
        .set_y_range(AutoOption::Fix(-1.1), AutoOption::Fix(2.5));

    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("qcd_mu{mu:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();

        let mut fout = io::BufWriter::new(
            fs::File::create(THIS_BASEDIR.join(dir).join(format!("qcd_mu{mu:.3}.out"))).unwrap(),
        );
        oms.iter().zip(spectr.iter()).for_each(|(om, s)| {
            writeln!(fout, "{om:.4}\t{s:.4}").unwrap();
        });
    }
}

fn main() {
    /* 0. Setup */
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    set_default_max_iter_integral(20);

    if !SHOW {
        fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
        fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();
    };

    let ommin = 0.02;
    let ommax = 2.0;
    let nom = 500;
    let dom = (ommax - ommin) / (nom as R);

    let oms: Vec<R> = (0..=nom).map(|i| ommin + (i as R) * dom).collect();

    /* I. Pure Yang-Mills */

    if !NO_YM {
        // IA. Fixed parameters
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("ym_fixed").as_path()).unwrap();
        }

        [
            0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
        ]
        .iter()
        .for_each(|&t| {
            plot_ym_t(&oms, t, MG, F0, "ym_fixed");
        });

        // IB. Parameters from the lattice
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("ym_lattice").as_path()).unwrap();
            fs::copy(
                THIS_BASEDIR.join("ym_fixed").join("ym_t0.000.png"),
                THIS_BASEDIR.join("ym_lattice").join("ym_t0.000.png"),
            )
            .unwrap();
        }

        [
            (0.121, 0.675, -0.83),
            (0.194, 0.725, -0.78),
            (0.260, 0.775, -0.58),
            (0.290, 0.725, -0.48),
            (0.366, 0.800, -0.38),
            (0.458, 0.900, -0.28),
        ]
        .iter()
        .for_each(|&(t, mg, f0)| {
            plot_ym_t(&oms, t, mg, f0, "ym_lattice");
        });
    }

    /* II. Full QCD, fixed parameters */

    if !NO_QCD {
        let fieldconfig = FieldConfig::new(NC, MG, vec![(2, 0.35), (1, 0.45)]);
        //let correctedfieldconfig = FieldConfig::new(NC, MG, vec![(2, 0.125), (1, 0.225)]);
        let correctedfieldconfig = fieldconfig.clone();

        // To avoid a spurious MOM-scheme divergence in the Mquark->0 limit
        // we must absorb ln(Mquark) into f0
        let f00 = F0
            + fieldconfig
                .quarks
                .iter()
                .map(|(nf, mq)| (*nf as R) * (fieldconfig.gluon / mq).ln())
                .sum::<R>()
                * 4.
                / (9. * (fieldconfig.nc as R));
        let f0c = f00
            + correctedfieldconfig
                .quarks
                .iter()
                .enumerate()
                .map(|(i, (nf, mqc))| {
                    let mq = fieldconfig.quarks[i].1;
                    (*nf as R) * (mq / mqc).ln()
                })
                .sum::<R>()
                * 4.
                / (9. * (fieldconfig.nc as R));

        // IIA. Zero density, as a function of temperature
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("qcd_zero_density_fixed").as_path()).unwrap();
        }

        [0., 0.05].iter().for_each(|&t| {
            plot_qcd_t(&oms, t, &fieldconfig, f00, "qcd_zero_density_fixed");
        });

        [
            0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
        ]
        .iter()
        .for_each(|&t| {
            plot_qcd_t(
                &oms,
                t,
                &correctedfieldconfig,
                f0c,
                "qcd_zero_density_fixed",
            );
        });

        // IIB. Zero temperature, as a function of chemical potential
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("qcd_zero_temperature_fixed").as_path()).unwrap();
        }

        [0., 0.1, 0.2, 0.3].iter().for_each(|&mu| {
            plot_qcd_mu(&oms, mu, &fieldconfig, f00, "qcd_zero_temperature_fixed");
        });

        [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].iter().for_each(|&mu| {
            plot_qcd_mu(
                &oms,
                mu,
                &correctedfieldconfig,
                f0c,
                "qcd_zero_temperature_fixed",
            );
        });

        // IIC. Zero density, as a function of temperature, Nf=2, lowest temperature params from lattice
        const MGQCD: R = 0.752;
        const F0QCD: R = -0.506;
        let fieldconfig = FieldConfig::new(NC, MGQCD, vec![(2, MQ)]);

        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("qcd_zero_density_fixed_nf2").as_path()).unwrap();
        }

        [
            0., 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7,
        ]
        .iter()
        .for_each(|&t| {
            plot_qcd_t(&oms, t, &fieldconfig, F0QCD, "qcd_zero_density_fixed_nf2");
        });
    }

    /* III. Full QCD, lattice parameters */

    if !NO_QCD {
        // IIIA. Zero density, as a function of temperature, Nf=2
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("qcd_zero_density_lattice").as_path()).unwrap();
        }

        // These data were obtained by renormalizing the zero-matsubara transverse propagator
        // at mu = 4 GeV, cutting the lattice data at 2 GeV and fixing mq = 400 MeV
        [
            (0., 0.7517777951910817, -0.5062071238487814), // This is a fake T = 0 value for reference
            (0.139, 0.7517777951910817, -0.5062071238487814),
            (0.154, 0.7639062870106536, -0.44084616449186326),
            (0.174, 0.7345150400347596, -0.3843769966635342),
            (0.199, 0.7306060280562273, -0.38161604237819363),
            (0.233, 0.746307345288752, -0.3585815488317115),
            (0.278, 0.7924693589714722, -0.32681035700549366),
        ]
        .iter()
        .for_each(|(t, mg, f0)| {
            let fieldconfig = FieldConfig::new(NC, *mg, vec![(2, MQ)]);
            plot_qcd_t(&oms, *t, &fieldconfig, *f0, "qcd_zero_density_lattice");
        });
    }
}
