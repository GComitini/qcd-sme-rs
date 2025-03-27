#[cfg(feature = "ftm_big")]
use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::qcd::thermal::gluon::{
    dressing_t_landau_w_field_config as qcd_dressing,
    dressing_t_zero_temp_landau_w_field_config as qcd_dressing_zero_temp,
};
use qcd_sme::qcd::FieldConfig;
use qcd_sme::types::Num;
use qcd_sme::types::{C, NCTYPE, R};
use qcd_sme::ym::gluon::dressing_landau as ym_dressing_zero_temp;
use qcd_sme::ym::thermal::gluon::dressing_t_landau as ym_dressing;
use tc_long_prop::BASEDIR;

use gnuplot::{AutoOption, AxesCommon, Figure, LabelOption, PlotOption};
use lazy_static::lazy_static;
use log::info;
use rayon::prelude::*;
use std::{
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
};

const NC: NCTYPE = 3;
const MG: R = 0.656;
const P0: R = 0.01;
const F0: R = -0.876;

const PIXX: u32 = 1600;
const PIXY: u32 = 1600;

const NO_QCD: bool = false;
const NO_YM: bool = false;
const SHOW: bool = false;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_poles");
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
                    { ym_dressing_zero_temp((om * om + P0 * P0) / (MG * MG), F0).abs() }
                    #[cfg(feature = "ftm_cut")]
                    { ym_dressing_zero_temp((om * om + P0 * P0) / (MG * MG), F0).im }
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
    info!(
        "Computed pure Yang-Mills dressing function at T = {t:.3} GeV. Pole is at ({omx:.4}) GeV."
    );

    #[cfg(not(feature = "ftm_cut"))]
    let (zmin, zmax) = (0., 100.);
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
                "T = {t:.3} GeV, \u{03c9} = ({z:.2}) GeV",
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

fn find_and_plot_qcd_t(
    oms: &[C],
    t: R,
    config: &FieldConfig,
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
                    { qcd_dressing_zero_temp(om, P0, 0., f0, config).abs() }
                    #[cfg(feature = "ftm_cut")]
                    { qcd_dressing_zero_temp(om, P0, 0., f0, config).im }
                };
                info!("Computed full QCD mu=0 dressing function at T = 0.000 GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1);
                res
            })
            .collect()
    } else {
        oms.par_iter().enumerate()
            .map(|(i, &om)| {
                let res = if om.re == 0. { R::NAN } else {
                    #[cfg(not(feature = "ftm_cut"))]
                    { qcd_dressing(om, P0, 1. / t, 0., f0, config).abs() }
                    #[cfg(feature = "ftm_cut")]
                    { qcd_dressing(om, P0, 1. / t, 0., f0, config).im }
                };
                info!("Computed full QCD mu=0 dressing function at T = {t:.3} GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1);
                res
            })
            .collect()
    };

    let omx = oms[max_om_idx(&mat)];
    info!("Computed full QCD mu=0 dressing function at T = {t:.3} GeV. Pole is at ({omx:.4}) GeV.");

    #[cfg(not(feature = "ftm_cut"))]
    let (zmin, zmax) = (0., 100.);
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
                "T = {t:.3} GeV, \u{03c9} = ({z:.2}) GeV",
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
                    .join(format!("qcd_t{t:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    omx
}

fn find_and_plot_qcd_mu(
    oms: &[C],
    mu: R,
    config: &FieldConfig,
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

    let mat: Vec<R> = oms
        .par_iter().enumerate()
        .map(|(i, &om)| {
            let res = if om.re == 0. { R::NAN } else {
                #[cfg(not(feature = "ftm_cut"))]
                { qcd_dressing_zero_temp(om, P0, mu, f0, config).abs() }
                #[cfg(feature = "ftm_cut")]
                { qcd_dressing_zero_temp(om, P0, mu, f0, config).im }
            };
            info!(
                "Computed full QCD T=0 dressing function at mu = {mu:.3} GeV for z = ({om:.4}) GeV ({}/{nact}): {res}", i+1
            );
            res
        })
        .collect();

    let omx = oms[max_om_idx(&mat)];
    info!(
        "Computed full QCD T=0 dressing function at mu = {mu:.3} GeV. Pole is at ({omx:.4}) GeV."
    );

    #[cfg(not(feature = "ftm_cut"))]
    let (zmin, zmax) = (0., 100.);
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
                "mu = {mu:.3} GeV, \u{03c9} = ({z:.2}) GeV",
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
                    .join(format!("qcd_mu{mu:.3}.png"))
                    .as_path(),
                PIXX,
                PIXY,
            )
            .unwrap();
    }

    omx
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    #[cfg(feature = "ftm_big")]
    set_default_max_iter_integral(20);

    if !SHOW {
        fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
        fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();
    };

    #[cfg(not(feature = "ftm_big"))]
    let ommax = 1.;
    #[cfg(feature = "ftm_big")]
    let ommax = 3.;

    #[cfg(not(feature = "ftm_twosided"))]
    let nom = 100;
    #[cfg(feature = "ftm_twosided")]
    let nom = 200;

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

    let (mut ym_fixed, mut ym_lattice, mut qcd_fixed_mu_zero, mut qcd_fixed_t_zero) =
        (vec![], vec![], vec![], vec![]);

    /* I. Pure Yang-Mills */

    if !NO_YM {
        // IA. Fixed parameters
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("ym_fixed").as_path()).unwrap();
        }

        [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
            .iter()
            .for_each(|&t| {
                let z = find_and_plot_ym_t(&oms, t, MG, F0, ommin, ommax, nom, "ym_fixed");
                ym_fixed.push((t, z.im, z.re))
            });
        let res0 = ym_fixed[0];

        let mut fpf = BufWriter::new(fs::File::create(THIS_BASEDIR.join("ym_fixed.out")).unwrap());
        ym_fixed
            .iter()
            .for_each(|(t, x, y)| writeln!(fpf, "{t:.4}\t{x:.4}\t{y:.4}").unwrap());

        // IB. Parameters from the lattice
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("ym_lattice").as_path()).unwrap();
            fs::copy(
                THIS_BASEDIR.join("ym_fixed").join("ym_t0.000.png"),
                THIS_BASEDIR.join("ym_lattice").join("ym_t0.000.png"),
            )
            .unwrap();
        }

        ym_lattice.push(res0);
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
            let z = find_and_plot_ym_t(&oms, t, mg, f0, ommin, ommax, nom, "ym_lattice");
            ym_lattice.push((t, z.im, z.re));
        });

        let mut fpf =
            BufWriter::new(fs::File::create(THIS_BASEDIR.join("ym_lattice.out")).unwrap());
        ym_lattice
            .iter()
            .for_each(|(t, x, y)| writeln!(fpf, "{t:.4}\t{x:.4}\t{y:.4}").unwrap());
    }

    /* II. Full QCD */

    if !NO_QCD {
        let fieldconfig = FieldConfig::new(NC, MG, vec![(2, 0.35), (1, 0.45)]);
        let correctedfieldconfig = FieldConfig::new(NC, MG, vec![(2, 0.125), (1, 0.225)]);

        // To avoid a spurious MOM-scheme divergence in the Mquark->0 limit
        // we must absorb ln(Mquark) into f0
        let f00 = F0
            + fieldconfig
                .quarks
                .iter()
                .map(|(nf, mq)| (*nf as R) * (MG / mq).ln())
                .sum::<R>()
                * 4.
                / (9. * (NC as R));
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

        // IIA. Zero density, as a function of temperature, fixed parameters
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("qcd_zero_density_fixed").as_path()).unwrap();
        }

        [0., 0.05].iter().for_each(|&t| {
            let z = find_and_plot_qcd_t(
                &oms,
                t,
                &fieldconfig,
                f00,
                ommin,
                ommax,
                nom,
                "qcd_zero_density_fixed",
            );
            qcd_fixed_mu_zero.push((t, z.im, z.re));
        });

        [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
            .iter()
            .for_each(|&t| {
                let z = find_and_plot_qcd_t(
                    &oms,
                    t,
                    &correctedfieldconfig,
                    f0c,
                    ommin,
                    ommax,
                    nom,
                    "qcd_zero_density_fixed",
                );
                qcd_fixed_mu_zero.push((t, z.im, z.re))
            });

        let mut fpf = BufWriter::new(
            fs::File::create(THIS_BASEDIR.join("qcd_zero_density_fixed.out")).unwrap(),
        );
        qcd_fixed_mu_zero
            .iter()
            .for_each(|(t, x, y)| writeln!(fpf, "{t:.4}\t{x:.4}\t{y:.4}").unwrap());

        // IIB. Zero temperature, as a function of chemical potential, fixed parameters
        if !SHOW {
            fs::create_dir_all(THIS_BASEDIR.join("qcd_zero_temperature_fixed").as_path()).unwrap();
        }

        [0., 0.1, 0.2, 0.3].iter().for_each(|&mu| {
            let z = find_and_plot_qcd_mu(
                &oms,
                mu,
                &fieldconfig,
                f00,
                ommin,
                ommax,
                nom,
                "qcd_zero_temperature_fixed",
            );
            qcd_fixed_t_zero.push((mu, z.im, z.re));
        });

        [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0].iter().for_each(|&mu| {
            let z = find_and_plot_qcd_mu(
                &oms,
                mu,
                &correctedfieldconfig,
                f0c,
                ommin,
                ommax,
                nom,
                "qcd_zero_temperature_fixed",
            );
            qcd_fixed_t_zero.push((mu, z.im, z.re))
        });

        let mut fpf = BufWriter::new(
            fs::File::create(THIS_BASEDIR.join("qcd_zero_temperature_fixed.out")).unwrap(),
        );
        qcd_fixed_t_zero
            .iter()
            .for_each(|(t, x, y)| writeln!(fpf, "{t:.4}\t{x:.4}\t{y:.4}").unwrap());
    }

    /* III. Pole trajectories */

    // IIIA. Ym fixed, ym lattice and zero-density qcd fixed, as a function of temperature
    if !(NO_YM && NO_QCD) {
        let mut figure = Figure::new();
        let ax2 = figure
            .axes2d()
            .set_x_label("Re(\u{03c9})", &[])
            .set_y_label("Im(\u{03c9})", &[])
            .set_x_range(AutoOption::Fix(0.5), AutoOption::Fix(1.))
            .set_y_range(AutoOption::Fix(0.3), AutoOption::Fix(0.8));

        if !NO_YM {
            ax2.lines_points(
                ym_fixed[0..ym_fixed.len() - 1]
                    .iter()
                    .map(|&(_, x, _)| x)
                    .collect::<Vec<R>>(),
                ym_fixed[0..ym_fixed.len() - 1]
                    .iter()
                    .map(|&(_, _, y)| y)
                    .collect::<Vec<R>>(),
                &[PlotOption::LineWidth(3.)],
            )
            .lines_points(
                ym_lattice.iter().map(|&(_, x, _)| x).collect::<Vec<R>>(),
                ym_lattice.iter().map(|&(_, _, y)| y).collect::<Vec<R>>(),
                &[PlotOption::LineWidth(3.)],
            );
        }

        if !NO_QCD {
            ax2.lines_points(
                qcd_fixed_mu_zero[0..ym_fixed.len() - 1]
                    .iter()
                    .map(|&(_, x, _)| x)
                    .collect::<Vec<R>>(),
                qcd_fixed_mu_zero[0..ym_fixed.len() - 1]
                    .iter()
                    .map(|&(_, _, y)| y)
                    .collect::<Vec<R>>(),
                &[PlotOption::LineWidth(3.)],
            );
        }

        if SHOW {
            figure.show().unwrap();
        } else {
            figure
                .save_to_png(THIS_BASEDIR.join("poles_t.png"), PIXX, PIXY)
                .unwrap();
        }
    }

    // IIIB. Zero-temperature qcd fixed
    if !NO_QCD {
        let mut figure = Figure::new();
        let ax2 = figure
            .axes2d()
            .set_x_label("Re(\u{03c9})", &[])
            .set_y_label("Im(\u{03c9})", &[])
            .set_x_range(AutoOption::Fix(0.5), AutoOption::Fix(1.))
            .set_y_range(AutoOption::Fix(0.3), AutoOption::Fix(0.8));

        ax2.lines_points(
            qcd_fixed_t_zero
                .iter()
                .map(|&(_, x, _)| x)
                .collect::<Vec<R>>(),
            qcd_fixed_t_zero
                .iter()
                .map(|&(_, _, y)| y)
                .collect::<Vec<R>>(),
            &[PlotOption::LineWidth(3.)],
        );

        if SHOW {
            figure.show().unwrap();
        } else {
            figure
                .save_to_png(THIS_BASEDIR.join("poles_mu.png"), PIXX, PIXY)
                .unwrap();
        }
    }
}
