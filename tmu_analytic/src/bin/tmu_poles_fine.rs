use argmin::core::observers::ObserverMode;
use argmin::core::{CostFunction, Executor, Gradient, IterState};
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin_observer_slog::SlogLogger;
use lazy_static::lazy_static;
use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::qcd::thermal::gluon::{
    dressing_t_inv_landau_w_field_config as qcd_dressing_inv,
    dressing_t_inv_zero_temp_landau_w_field_config as qcd_dressing_inv_zero_temp,
};
use qcd_sme::qcd::FieldConfig;
use qcd_sme::types::{C, NCTYPE, R};
use qcd_sme::ym::gluon::dressing_inv_landau as ym_dressing_inv_zero_temp;
use qcd_sme::ym::thermal::gluon::dressing_t_inv_landau as ym_dressing_inv;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use tmu_analytic::BASEDIR;

const EPS_GRAD: R = 1E-6;
const MAX_SD_ITER: u64 = 5000;
const P0: R = 0.01;
const TARGET_ZERO: R = 1E-12;

const NC: NCTYPE = 3;
const MG: R = 0.656;
const F0: R = -0.876;
const MQ: R = 0.4;

const NO_YM: bool = false;
const NO_QCD: bool = false;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_poles_fine");
}

pub struct YMZeroFinder {
    pub mg: R,
    pub f0: R,
    pub t: R,
}

fn cost_ym(z: C, t: R, m: R, f0: R) -> R {
    let r = if t == 0. {
        let s = (z * z + P0 * P0) / (m * m);
        ym_dressing_inv_zero_temp(s, f0).norm()
    } else {
        ym_dressing_inv(z, P0, m, 1. / t, f0).norm()
    };
    r * r
}

impl CostFunction for YMZeroFinder {
    type Param = Vec<R>;

    type Output = R;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let z = C::new(param[0], param[1]);
        Ok(cost_ym(z, self.t, self.mg, self.f0))
    }
}

impl Gradient for YMZeroFinder {
    type Param = Vec<R>;

    type Gradient = Vec<R>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, argmin::core::Error> {
        let (x, y) = (param[0], param[1]);
        Ok(vec![
            (cost_ym(x + EPS_GRAD + C::I * y, self.t, self.mg, self.f0)
                - cost_ym(x - EPS_GRAD + C::I * y, self.t, self.mg, self.f0))
                / (2. * EPS_GRAD),
            (cost_ym(x + C::I * (y + EPS_GRAD), self.t, self.mg, self.f0)
                - cost_ym(x + C::I * (y - EPS_GRAD), self.t, self.mg, self.f0))
                / (2. * EPS_GRAD),
        ])
    }
}

#[allow(clippy::type_complexity)]
impl YMZeroFinder {
    pub fn executor(
        self,
        init_param: Vec<R>,
        observe: bool,
    ) -> Executor<
        YMZeroFinder,
        SteepestDescent<MoreThuenteLineSearch<Vec<R>, Vec<R>, R>>,
        IterState<Vec<R>, Vec<R>, (), (), (), R>,
    > {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);
        let mut exec = Executor::new(self, solver).configure(|state| {
            state
                .param(init_param)
                .max_iters(MAX_SD_ITER)
                .target_cost(TARGET_ZERO)
        });
        if observe {
            exec = exec.add_observer(SlogLogger::term(), ObserverMode::Always)
        }
        exec
    }
}

pub struct QCDZeroFinder {
    pub config: FieldConfig,
    pub f0: R,
    pub t: R,
}

fn cost_qcd(z: C, t: R, f0: R, config: &FieldConfig) -> R {
    let r = if t == 0. {
        qcd_dressing_inv_zero_temp(z, P0, 0., f0, config).norm()
    } else {
        qcd_dressing_inv(z, P0, 1. / t, 0., f0, config).norm()
    };
    r * r
}

impl CostFunction for QCDZeroFinder {
    type Param = Vec<R>;

    type Output = R;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let z = C::new(param[0], param[1]);
        Ok(cost_qcd(z, self.t, self.f0, &self.config))
    }
}

impl Gradient for QCDZeroFinder {
    type Param = Vec<R>;

    type Gradient = Vec<R>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, argmin::core::Error> {
        let (x, y) = (param[0], param[1]);
        Ok(vec![
            (cost_qcd(x + EPS_GRAD + C::I * y, self.t, self.f0, &self.config)
                - cost_qcd(x - EPS_GRAD + C::I * y, self.t, self.f0, &self.config))
                / (2. * EPS_GRAD),
            (cost_qcd(x + C::I * (y + EPS_GRAD), self.t, self.f0, &self.config)
                - cost_qcd(x + C::I * (y - EPS_GRAD), self.t, self.f0, &self.config))
                / (2. * EPS_GRAD),
        ])
    }
}

#[allow(clippy::type_complexity)]
impl QCDZeroFinder {
    pub fn executor(
        self,
        init_param: Vec<R>,
        observe: bool,
    ) -> Executor<
        QCDZeroFinder,
        SteepestDescent<MoreThuenteLineSearch<Vec<R>, Vec<R>, R>>,
        IterState<Vec<R>, Vec<R>, (), (), (), R>,
    > {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);
        let mut exec = Executor::new(self, solver).configure(|state| {
            state
                .param(init_param)
                .max_iters(MAX_SD_ITER)
                .target_cost(TARGET_ZERO)
        });
        if observe {
            exec = exec.add_observer(SlogLogger::term(), ObserverMode::Always)
        }
        exec
    }
}

fn main() {
    set_default_max_iter_integral(20);

    fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
    fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();

    /* I. Pure Yang-Mills */
    if !NO_YM {
        let (mut ym_fixed, mut ym_lattice) = (vec![], vec![]);

        // IA. Fixed parameters
        let mut prev_zero = vec![1., 1.];

        let temps = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
        temps.iter().for_each(|&t| {
            println!("Searching YM zero at fixed parameters for T = {t:.3}");
            let ymzf = YMZeroFinder { mg: MG, f0: F0, t };
            let res = ymzf.executor(prev_zero.clone(), true).run().unwrap();
            println!("{res}");
            let zero = res.state.best_param.unwrap();
            ym_fixed.push((vec![zero[1], zero[0]], res.state.best_cost));
            prev_zero = zero;
        });

        let mut fout = BufWriter::new(File::create(THIS_BASEDIR.join("ym_fixed.out")).unwrap());
        temps.iter().zip(ym_fixed.iter()).for_each(|(t, pole)| {
            writeln!(
                fout,
                "{t:.3}\t{:.5}\t{:.5}\t{:.24}",
                pole.0[0], pole.0[1], pole.1
            )
            .unwrap()
        });

        // IB. Lattice parameters
        let mut prev_zero = vec![1., 1.];

        let params = [
            (0.121, 0.675, -0.83),
            (0.194, 0.725, -0.78),
            (0.260, 0.775, -0.58),
            (0.290, 0.725, -0.48),
            (0.366, 0.800, -0.38),
            (0.458, 0.900, -0.28),
        ];
        params.iter().for_each(|&(t, mg, f0)| {
            println!("Searching YM zero at lattice parameters for T = {t:.3}");
            let ymzf = YMZeroFinder { mg, f0, t };
            let res = ymzf.executor(prev_zero.clone(), true).run().unwrap();
            println!("{res}");
            let zero = res.state.best_param.unwrap();
            ym_lattice.push((vec![zero[1], zero[0]], res.state.best_cost));
            prev_zero = zero;
        });

        let mut fout = BufWriter::new(File::create(THIS_BASEDIR.join("ym_lattice.out")).unwrap());
        params
            .iter()
            .zip(ym_lattice.iter())
            .for_each(|((t, _, _), pole)| {
                writeln!(
                    fout,
                    "{t:.3}\t{:.5}\t{:.5}\t{:.24}",
                    pole.0[0], pole.0[1], pole.1
                )
                .unwrap()
            });
    }

    /* II. QCD */
    if !NO_QCD {
        let (mut qcd_fixed, mut qcd_fixed_nf2, mut qcd_lattice) = (vec![], vec![], vec![]);

        // IIA. Fixed parameters
        let mut prev_zero = vec![-1., 1.];

        let fieldconfig = FieldConfig::new(NC, MG, vec![(2, 0.35), (1, 0.45)]);
        let f0q = F0
            + fieldconfig
                .quarks
                .iter()
                .map(|(nf, mq)| (*nf as R) * (fieldconfig.gluon / mq).ln())
                .sum::<R>()
                * 4.
                / (9. * (fieldconfig.nc as R));
        let temps = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
        temps.iter().for_each(|&t| {
            println!("Searching QCD zero at fixed parameters for T = {t:.3}");
            let qcdzf = QCDZeroFinder {
                config: fieldconfig.clone(),
                f0: f0q,
                t,
            };
            let res = qcdzf.executor(prev_zero.clone(), true).run().unwrap();
            println!("{res}");
            let zero = res.state.best_param.unwrap();
            qcd_fixed.push((vec![zero[1], zero[0]], res.state.best_cost));
            prev_zero = zero;
        });

        let mut fout = BufWriter::new(File::create(THIS_BASEDIR.join("qcd_fixed.out")).unwrap());
        temps.iter().zip(qcd_fixed.iter()).for_each(|(t, pole)| {
            writeln!(
                fout,
                "{t:.3}\t{:.5}\t{:.5}\t{:.24}",
                pole.0[0], pole.0[1], pole.1
            )
            .unwrap()
        });

        // IIB. Nf = 2, lowest temperature fixed params from lattice
        let mut prev_zero = vec![1., 1.];

        const MGQCD: R = 0.752;
        const F0QCD: R = -0.506;
        let fieldconfig = FieldConfig::new(NC, MGQCD, vec![(2, MQ)]);
        let temps = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
        temps.iter().for_each(|&t| {
            println!("Searching Nf = 2 QCD zero at fixed parameters for T = {t:.3}");
            let qcdzf = QCDZeroFinder {
                config: fieldconfig.clone(),
                f0: F0QCD,
                t,
            };
            let res = qcdzf.executor(prev_zero.clone(), true).run().unwrap();
            println!("{res}");
            let zero = res.state.best_param.unwrap();
            qcd_fixed_nf2.push((vec![zero[1], zero[0]], res.state.best_cost));
            prev_zero = zero;
        });

        let mut fout =
            BufWriter::new(File::create(THIS_BASEDIR.join("qcd_fixed_nf2.out")).unwrap());
        temps
            .iter()
            .zip(qcd_fixed_nf2.iter())
            .for_each(|(t, pole)| {
                writeln!(
                    fout,
                    "{t:.3}\t{:.5}\t{:.5}\t{:.24}",
                    pole.0[0], pole.0[1], pole.1
                )
                .unwrap()
            });

        // IIC. Lattice parameters
        let mut prev_zero = vec![1., 1.];

        let params = [
            (0.139, 0.7517777951910817, -0.5062071238487814),
            (0.154, 0.7639062870106536, -0.44084616449186326),
            (0.174, 0.7345150400347596, -0.3843769966635342),
            (0.199, 0.7306060280562273, -0.38161604237819363),
            (0.233, 0.746307345288752, -0.3585815488317115),
            (0.278, 0.7924693589714722, -0.32681035700549366),
        ];
        params.iter().for_each(|&(t, mg, f0)| {
            println!("Searching QCD zero at lattice parameters for T = {t:.3}");
            let fieldconfig = FieldConfig::new(NC, mg, vec![(2, MQ)]);
            let qcdzf = QCDZeroFinder {
                config: fieldconfig,
                f0,
                t,
            };
            let res = qcdzf.executor(prev_zero.clone(), true).run().unwrap();
            println!("{res}");
            let zero = res.state.best_param.unwrap();
            qcd_lattice.push((vec![zero[1], zero[0]], res.state.best_cost));
            prev_zero = zero;
        });

        let mut fout = BufWriter::new(File::create(THIS_BASEDIR.join("qcd_lattice.out")).unwrap());
        params
            .iter()
            .zip(qcd_lattice.iter())
            .for_each(|((t, _, _), pole)| {
                writeln!(
                    fout,
                    "{t:.3}\t{:.5}\t{:.5}\t{:.24}",
                    pole.0[0], pole.0[1], pole.1
                )
                .unwrap()
            });
    }
}
