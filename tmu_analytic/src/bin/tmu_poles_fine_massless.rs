use argmin::core::observers::ObserverMode;
use argmin::core::{CostFunction, Executor, Gradient, IterState};
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin_observer_slog::SlogLogger;
use lazy_static::lazy_static;
use qcd_sme::consts::set_default_max_iter_integral;
use qcd_sme::types::{C, R};
use qcd_sme::ym::gluon::dressing_inv_landau as ym_dressing_inv_zero_temp;
use qcd_sme::ym::thermal::gluon::dressing_t_inv_landau as ym_dressing_inv;
use std::fs;
use std::path::PathBuf;
use tmu_analytic::BASEDIR;

const EPS_GRAD: R = 1E-6;
const MAX_SD_ITER: u64 = 5000;
const TARGET_ZERO: R = 1E-12;

const MG: R = 0.005;
const P0: R = 0.01;
const F0: R = 25.;
const T: R = 0.1;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("fintmu_poles_fine_massless");
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

fn main() {
    set_default_max_iter_integral(20);

    fs::remove_dir_all(THIS_BASEDIR.as_path()).ok();
    fs::create_dir_all(THIS_BASEDIR.as_path()).unwrap();

    let prev_zero = vec![0.013865064, 0.0417031146];

    println!("Searching YM zero at fixed parameters for T = {T:.3}, F0 = {F0}");
    let ymzf = YMZeroFinder {
        mg: MG,
        f0: F0,
        t: T,
    };
    let res = ymzf.executor(prev_zero.clone(), true).run().unwrap();
    println!("{res}");
    let zero = res.state.best_param.unwrap();
    println!("");
    println!("POLE: {}, {}, {}", zero[1], zero[0], res.state.best_cost);
    // POLE: 0.04170311460987214, 0.013865064016310102, 0.0000000000000009183932171457694
}
