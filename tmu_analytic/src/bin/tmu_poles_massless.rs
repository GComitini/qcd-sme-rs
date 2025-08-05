#[cfg(feature = "ftm_big")]
use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::types::Num;
use qcd_sme::types::{C, R};
use qcd_sme::ym::gluon::dressing_landau as ym_dressing_zero_temp;
use qcd_sme::ym::thermal::gluon::dressing_t_landau as ym_dressing;
use tmu_analytic::BASEDIR;

use gnuplot::{AutoOption, AxesCommon, Figure, LabelOption};
use lazy_static::lazy_static;
use log::info;
use rayon::prelude::*;
use std::{
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
};

const MG: R = 0.005;
const P0: R = 0.01;
const F0: R = 25.;
#[cfg(feature = "ftm_proptest")]
const PREN: R = 0.5;
#[cfg(feature = "ftm_proptest")]
const OMEPS: R = 1E-3;

const PIXX: u32 = 1600;
const PIXY: u32 = 1600;

const SHOW: bool = false;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_poles_massless");
}

fn max_om_idx(mat: &[R]) -> usize {
    let (mut mx, mut imx) = (0., 0);
    for (i, &m) in mat.iter().enumerate() {
        if m > mx {
            mx = m;
            imx = i;
        }
    }
    imx
}

#[allow(clippy::too_many_arguments)]
fn find_and_plot_ym_t(
    oms: &[C],
    t: R,
    m: R,
    f0: R,
    ommin: R,
    ommax: R,
    nom: usize,
    dir: &str,
) -> C {
    #[cfg(not(feature = "ftm_twosided"))]
    let nrc = nom;
    #[cfg(feature = "ftm_twosided")]
    let nrc = nom + 1;

    let nact = nrc * nrc;

    let mat: Vec<R> = if t == 0. {
        oms.par_iter().enumerate()
            .map(|(i, &om)| {
                let res = if om.re == 0. { R::NAN } else {
                    #[cfg(not(feature = "ftm_cut"))]
                    { ym_dressing_zero_temp((om * om + P0 * P0) / (m*m), f0).abs() }
                    #[cfg(feature = "ftm_cut")]
                    { ym_dressing_zero_temp((om * om + P0 * P0) / (m*m), f0).im }
                };
                info!(
                    "Computed pure Yang-Mills dressing function at T = 0.000 GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1
                );
                res
            })
            .collect()
    } else {
        oms.par_iter().enumerate()
            .map(|(i,&om)| {
                let res = if om.re == 0. { R::NAN } else {
                    #[cfg(not(feature = "ftm_cut"))]
                    { ym_dressing(om, P0, m, 1. / t, f0).abs() }
                    #[cfg(feature = "ftm_cut")]
                    { ym_dressing(om, P0, m, 1. / t, f0).im }
                };
                info!(
                    "Computed pure Yang-Mills dressing function at T = {t:.3} GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1
                );
                res
            })
            .collect()
    };

    let omx = oms[max_om_idx(&mat)];
    let omx = C::new(omx.re.abs(), omx.im.abs());
    info!(
        "Computed pure Yang-Mills dressing function at T = {t:.3} GeV. Pole is at ({omx:.4}) GeV."
    );

    #[cfg(not(feature = "ftm_cut"))]
    let (zmin, zmax) = (0., 8.);
    #[cfg(feature = "ftm_cut")]
    let (zmin, zmax) = (-2., 2.);
    let mut figure = Figure::new();
    figure
        .axes3d()
        .surface(mat, nrc, nrc, Some((ommin, ommin, ommax, ommax)), &[])
        .set_x_range(AutoOption::Fix(ommin), AutoOption::Fix(ommax))
        .set_y_range(AutoOption::Fix(ommin), AutoOption::Fix(ommax))
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
                    .join(format!("ym_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    omx
}

// If normalize is true, input momenta are treated as adimensionalized by MG and
// the output propagator is adimensionalized by MG. Note that MG is a fixed scale,
// not the (possibly temperature-dependent) mass parameter at the given temperature
#[cfg(feature = "ftm_proptest")]
fn ym_t_proptest(momenta: &[R], t: R, m: R, f0: R, normalize: bool, dir: &str) {
    let prop: Vec<R> = if t == 0. {
        let mut z = 1. / ym_dressing_zero_temp((OMEPS * OMEPS + PREN * PREN) / (m * m), f0).re();
        if normalize {
            z *= MG * MG;
        }
        momenta
            .par_iter()
            .map(|&p| {
                let p2 = if normalize {p*p*MG*MG} else {p*p};
                let res = z * ym_dressing_zero_temp((OMEPS * OMEPS + p2) / (m * m), f0).re()
                    / (OMEPS * OMEPS + p2);
                info!(
                    "PROPTEST: Computed pure Yang-Mills propagator at T = 0.000 GeV for p = {p:.4} GeV: {res}"
                );
                res
            })
            .collect()
    } else {
        let mut z = 1. / ym_dressing(OMEPS, PREN, m, 1. / t, f0).re;
        if normalize {
            z *= MG * MG;
        };
        momenta.par_iter()
            .map(|&p| {
                let p = if normalize {p*MG} else {p};
                let res = z * ym_dressing(OMEPS, p, m, 1. / t, f0).re/(OMEPS * OMEPS + p * p);
                info!(
                    "PROPTEST: Computed pure Yang-Mills propagator at T = {t:.3} GeV for p = {p:.4} GeV: {res}"
                );
                res
            })
            .collect()
    };

    let mut figure = Figure::new();
    let zmax = 5500.;
    figure
        .axes2d()
        .lines(momenta, prop, &[])
        .set_x_range(
            AutoOption::Fix(*momenta.first().unwrap()),
            AutoOption::Fix(*momenta.last().unwrap()),
        )
        .set_y_range(AutoOption::Fix(0.), AutoOption::Fix(zmax));
    if SHOW {
        figure.show().unwrap();
    } else {
        figure
            .save_to_png(
                THIS_BASEDIR
                    .join(dir)
                    .join(format!("proptest_ym_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }
}

fn main() {
    /* 0. Setup */
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    #[cfg(feature = "ftm_big")]
    set_default_max_iter_integral(20);

    if !SHOW {
        fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
        fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();
    };

    #[cfg(not(feature = "ftm_big"))]
    let ommax = 0.05;
    #[cfg(feature = "ftm_big")]
    let ommax = 3.;

    #[cfg(not(feature = "ftm_twosided"))]
    let nom = 500;
    #[cfg(feature = "ftm_twosided")]
    let nom = 500;

    #[cfg(not(feature = "ftm_twosided"))]
    let nom_actual = nom * nom;
    #[cfg(feature = "ftm_twosided")]
    let nom_actual = nom * nom + 2 * nom + 1;

    #[cfg(not(feature = "ftm_twosided"))]
    let dom = ommax / (nom as R);
    #[cfg(feature = "ftm_twosided")]
    let dom = 2. * ommax / (nom as R);

    #[cfg(not(feature = "ftm_twosided"))]
    let ommin = dom;
    #[cfg(feature = "ftm_twosided")]
    let ommin = -ommax;

    let mut oms = Vec::with_capacity(nom_actual);

    #[cfg(not(feature = "ftm_twosided"))]
    for i in 1..=nom {
        for j in 1..=nom {
            oms.push(C::new((j as R) * dom, (i as R) * dom));
        }
    }

    #[cfg(feature = "ftm_twosided")]
    for i in 0..=nom {
        for j in 0..=nom {
            oms.push(C::new((j as R) * dom - ommax, (i as R) * dom - ommax));
        }
    }

    let mut ym_fixed = vec![];

    #[cfg(feature = "ftm_proptest")]
    let ym_momenta: Vec<R> = {
        let (pmin, pmax, np) = (0.01, 0.2, 200);
        let dp = (pmax - pmin) / (np as R);
        (0..=np).map(|i| pmin + (i as R) * dp).collect()
    };

    if !SHOW {
        fs::create_dir_all(THIS_BASEDIR.join("ym_fixed").as_path()).unwrap();
    }

    let t = 0.1;
    let z = find_and_plot_ym_t(&oms, t, MG, F0, ommin, ommax, nom, "ym_fixed");
    ym_fixed.push((t, z.im, z.re));
    #[cfg(feature = "ftm_proptest")]
    ym_t_proptest(&ym_momenta, t, MG, F0, false, "ym_fixed");

    let mut fpf = BufWriter::new(fs::File::create(THIS_BASEDIR.join("ym_fixed.out")).unwrap());
    ym_fixed
        .iter()
        .for_each(|(t, x, y)| writeln!(fpf, "{t:.4}\t{x:.4}\t{y:.4}").unwrap());

    info!("POLES: {ym_fixed:?}");
}
