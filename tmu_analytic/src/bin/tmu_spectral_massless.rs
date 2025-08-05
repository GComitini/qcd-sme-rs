use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::types::{C, R};
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

const MG: R = 0.005;
const P0: R = 0.01;
const F0: R = 25.;
const PREN: R = 0.5;
const OMEPS: R = 1E-3;

const PIXX: u32 = 1600;
const PIXY: u32 = 1600;

const SHOW: bool = false;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_spectral_massless");
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
        .set_y_range(AutoOption::Fix(-300.), AutoOption::Fix(100.));

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

fn main() {
    /* 0. Setup */
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    set_default_max_iter_integral(20);

    if !SHOW {
        fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
        fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();
    };

    let ommin = 0.02;
    let ommax = 0.2;
    let nom = 200;
    let dom = (ommax - ommin) / (nom as R);

    let oms: Vec<R> = (0..=nom).map(|i| ommin + (i as R) * dom).collect();

    /* I. Pure Yang-Mills */
    if !SHOW {
        fs::create_dir_all(THIS_BASEDIR.join("ym_fixed").as_path()).unwrap();
    }

    [0., 0.05, 0.1, 0.15, 0.2, 0.25].iter().for_each(|&t| {
        plot_ym_t(&oms, t, MG, F0, "ym_fixed");
    });
}
