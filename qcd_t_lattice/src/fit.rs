use argmin::core::observers::ObserverMode;
use argmin::core::{CostFunction, Executor, Gradient, IterState};
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin_observer_slog::SlogLogger;
use qcd_sme::R;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

const EPS_GRAD: R = 1E-6;
const MAX_SD_ITER: u64 = 5000;

// Returns the reduced chi square computed for parameters (a, b) given the data
// [(p, d, e)] (independent variable, dependent variable, error) and the function
// ff(p, a, b)
fn cost2<F: Fn(R, R, R) -> R + Sync>(a: R, b: R, data: &[(R, R, R)], ff: &F) -> R {
    let diff2s: Vec<R> = data
        .par_iter()
        .map(|(p, d, e)| {
            let val = ff(*p, a, b);
            if val.is_nan() {
                // Numerical integrals can evaluate to NaN. If they do, assign
                // them a difference of 2e
                return 4.;
            };
            let diff = (val - d) / e;
            diff * diff
        })
        .collect();
    let mut chi = 0.;
    diff2s.iter().for_each(|df| chi += df);
    chi / ((data.len() - 2) as R)
}

pub struct ChiSquare2Params<'a, F: Fn(R, R, R) -> R> {
    pub ff: &'a F,
    pub data: &'a [(R, R, R)],
}

impl<F: Fn(R, R, R) -> R + Sync> CostFunction for ChiSquare2Params<'_, F> {
    type Param = Vec<R>;

    type Output = R;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok(cost2(param[0], param[1], self.data, self.ff))
    }
}

impl<F: Fn(R, R, R) -> R + Sync> Gradient for ChiSquare2Params<'_, F> {
    type Param = Vec<R>;

    type Gradient = Vec<R>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, argmin::core::Error> {
        let (a, b) = (param[0], param[1]);
        Ok(vec![
            (cost2(a + EPS_GRAD, b, self.data, self.ff)
                - cost2(a - EPS_GRAD, b, self.data, self.ff))
                / (2. * EPS_GRAD),
            (cost2(a, b + EPS_GRAD, self.data, self.ff)
                - cost2(a, b - EPS_GRAD, self.data, self.ff))
                / (2. * EPS_GRAD),
        ])
    }
}

#[allow(clippy::type_complexity)]
impl<'a, F: Fn(R, R, R) -> R + Sync> ChiSquare2Params<'a, F> {
    pub fn executor(
        self,
        init_param: Vec<R>,
        observe: bool,
    ) -> Executor<
        ChiSquare2Params<'a, F>,
        SteepestDescent<MoreThuenteLineSearch<Vec<R>, Vec<R>, R>>,
        IterState<Vec<R>, Vec<R>, (), (), (), R>,
    > {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);
        let mut exec =
            Executor::new(self, solver).configure(|state| state.param(init_param).max_iters(MAX_SD_ITER));
        if observe {
            exec = exec.add_observer(SlogLogger::term(), ObserverMode::Always)
        }
        exec
    }
}

pub fn paramgradient2<F: Fn(R, R, R) -> R + Sync>(f: F, p: R, a: R, b: R) -> Vec<R> {
    vec![
        (f(p, a + EPS_GRAD, b) - f(p, a - EPS_GRAD, b)) / (2. * EPS_GRAD),
        (f(p, a, b + EPS_GRAD) - f(p, a, b - EPS_GRAD)) / (2. * EPS_GRAD),
    ]
}

// Returns the reduced chi square computed for parameters (a, b, c) given the data
// [(p, d, e)] (independent variable, dependent variable, error) and the function
// ff(p, a, b, c)
fn cost3<F: Fn(R, R, R, R) -> R + Sync>(a: R, b: R, c: R, data: &[(R, R, R)], ff: &F) -> R {
    let diff2s: Vec<R> = data
        .par_iter()
        .map(|(p, d, e)| {
            let val = ff(*p, a, b, c);
            if val.is_nan() {
                // Numerical integrals can evaluate to NaN. If they do, assign
                // them a difference of 2e
                return 4.;
            };
            let diff = (val - d) / e;
            diff * diff
        })
        .collect();
    let mut chi = 0.;
    diff2s.iter().for_each(|df| chi += df);
    chi / ((data.len() - 3) as R)
}

pub struct ChiSquare3Params<'a, F: Fn(R, R, R, R) -> R + Sync> {
    pub ff: &'a F,
    pub data: &'a [(R, R, R)],
}

impl<F: Fn(R, R, R, R) -> R + Sync> CostFunction for ChiSquare3Params<'_, F> {
    type Param = Vec<R>;

    type Output = R;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok(cost3(param[0], param[1], param[2], self.data, self.ff))
    }
}

impl<F: Fn(R, R, R, R) -> R + Sync> Gradient for ChiSquare3Params<'_, F> {
    type Param = Vec<R>;

    type Gradient = Vec<R>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, argmin::core::Error> {
        let (a, b, c) = (param[0], param[1], param[2]);
        Ok(vec![
            (cost3(a + EPS_GRAD, b, c, self.data, self.ff)
                - cost3(a - EPS_GRAD, b, c, self.data, self.ff))
                / (2. * EPS_GRAD),
            (cost3(a, b + EPS_GRAD, c, self.data, self.ff)
                - cost3(a, b - EPS_GRAD, c, self.data, self.ff))
                / (2. * EPS_GRAD),
            (cost3(a, b, c + EPS_GRAD, self.data, self.ff)
                - cost3(a, b, c - EPS_GRAD, self.data, self.ff))
                / (2. * EPS_GRAD),
        ])
    }
}

#[allow(clippy::type_complexity)]
impl<'a, F: Fn(R, R, R, R) -> R + Sync> ChiSquare3Params<'a, F> {
    pub fn executor(
        self,
        init_param: Vec<R>,
        observe: bool,
    ) -> Executor<
        ChiSquare3Params<'a, F>,
        SteepestDescent<MoreThuenteLineSearch<Vec<R>, Vec<R>, R>>,
        IterState<Vec<R>, Vec<R>, (), (), (), R>,
    > {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);
        let mut exec =
            Executor::new(self, solver).configure(|state| state.param(init_param).max_iters(MAX_SD_ITER));
        if observe {
            exec = exec.add_observer(SlogLogger::term(), ObserverMode::Always)
        }
        exec
    }
}

pub fn paramgradient3<F: Fn(R, R, R, R) -> R + Sync>(f: F, p: R, a: R, b: R, c: R) -> Vec<R> {
    vec![
        (f(p, a + EPS_GRAD, b, c) - f(p, a - EPS_GRAD, b, c)) / (2. * EPS_GRAD),
        (f(p, a, b + EPS_GRAD, c) - f(p, a, b - EPS_GRAD, c)) / (2. * EPS_GRAD),
        (f(p, a, b, c + EPS_GRAD) - f(p, a, b, c - EPS_GRAD)) / (2. * EPS_GRAD),
    ]
}

// Returns the reduced chi square computed for parameters (a, b, c, d) given the data
// [(p, dd, e)] (independent variable, dependent variable, error) and the function
// ff(p, a, b, c, d)
fn cost4<F: Fn(R, R, R, R, R) -> R + Sync>(
    a: R,
    b: R,
    c: R,
    d: R,
    data: &[(R, R, R)],
    ff: &F,
) -> R {
    let diff4s: Vec<R> = data
        .par_iter()
        .map(|(p, dd, e)| {
            let val = ff(*p, a, b, c, d);
            if val.is_nan() {
                // Numerical integrals can evaluate to NaN. If they do, assign
                // them a difference of 2e
                return 4.;
            };
            let diff = (val - dd) / e;
            diff * diff
        })
        .collect();
    let mut chi = 0.;
    diff4s.iter().for_each(|df| chi += df);
    chi / ((data.len() - 4) as R)
}

pub struct ChiSquare4Params<'a, F: Fn(R, R, R, R, R) -> R + Sync> {
    pub ff: &'a F,
    pub data: &'a [(R, R, R)],
}

impl<F: Fn(R, R, R, R, R) -> R + Sync> CostFunction for ChiSquare4Params<'_, F> {
    type Param = Vec<R>;

    type Output = R;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok(cost4(
            param[0], param[1], param[2], param[3], self.data, self.ff,
        ))
    }
}

impl<F: Fn(R, R, R, R, R) -> R + Sync> Gradient for ChiSquare4Params<'_, F> {
    type Param = Vec<R>;

    type Gradient = Vec<R>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, argmin::core::Error> {
        let (a, b, c, d) = (param[0], param[1], param[2], param[3]);
        Ok(vec![
            (cost4(a + EPS_GRAD, b, c, d, self.data, self.ff)
                - cost4(a - EPS_GRAD, b, c, d, self.data, self.ff))
                / (2. * EPS_GRAD),
            (cost4(a, b + EPS_GRAD, c, d, self.data, self.ff)
                - cost4(a, b - EPS_GRAD, c, d, self.data, self.ff))
                / (2. * EPS_GRAD),
            (cost4(a, b, c + EPS_GRAD, d, self.data, self.ff)
                - cost4(a, b, c - EPS_GRAD, d, self.data, self.ff))
                / (2. * EPS_GRAD),
            (cost4(a, b, c, d + EPS_GRAD, self.data, self.ff)
                - cost4(a, b, c, d - EPS_GRAD, self.data, self.ff))
                / (2. * EPS_GRAD),
        ])
    }
}

#[allow(clippy::type_complexity)]
impl<'a, F: Fn(R, R, R, R, R) -> R + Sync> ChiSquare4Params<'a, F> {
    pub fn executor(
        self,
        init_param: Vec<R>,
        observe: bool,
    ) -> Executor<
        ChiSquare4Params<'a, F>,
        SteepestDescent<MoreThuenteLineSearch<Vec<R>, Vec<R>, R>>,
        IterState<Vec<R>, Vec<R>, (), (), (), R>,
    > {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = SteepestDescent::new(linesearch);
        let mut exec =
            Executor::new(self, solver).configure(|state| state.param(init_param).max_iters(MAX_SD_ITER));
        if observe {
            exec = exec.add_observer(SlogLogger::term(), ObserverMode::Always)
        }
        exec
    }
}

pub fn gradient4<F: Fn(R, R, R, R, R) -> R + Sync>(f: F, p: R, a: R, b: R, c: R, d: R) -> Vec<R> {
    vec![
        (f(p, a + EPS_GRAD, b, c, d) - f(p, a - EPS_GRAD, b, c, d)) / (2. * EPS_GRAD),
        (f(p, a, b + EPS_GRAD, c, d) - f(p, a, b - EPS_GRAD, c, d)) / (2. * EPS_GRAD),
        (f(p, a, b, c + EPS_GRAD, d) - f(p, a, b, c - EPS_GRAD, d)) / (2. * EPS_GRAD),
        (f(p, a, b, c, d + EPS_GRAD) - f(p, a, b, c, d - EPS_GRAD)) / (2. * EPS_GRAD),
    ]
}

pub fn readfile<P: AsRef<std::path::Path>>(p: P) -> Vec<(R, R, R)> {
    let mut fs = BufReader::new(File::open(p).unwrap());
    let mut buf = String::new();

    let mut data = Vec::new();

    while let Ok(l) = fs.read_line(&mut buf) {
        if l == 0 {
            break;
        }
        let linebuf: Vec<R> = buf
            .split_whitespace()
            .map(|a| a.parse::<R>().unwrap())
            .collect();
        data.push((linebuf[0], linebuf[1], linebuf[2]));
        buf.clear();
    }

    data
}

pub fn param_errors(data: &[(R, R, R)], grad: &[Vec<R>], chi2: R) -> Vec<R> {
    assert_eq!(data.len(), grad.len());
    let dim = grad.first().unwrap().len();
    let mut mat = nalgebra::DMatrix::<R>::zeros(dim, dim);
    data.iter().zip(grad.iter()).for_each(|((_, _, e), g)| {
        let gg: Vec<R> = g.to_vec();
        let gg = nalgebra::DVector::from_vec(gg);
        mat += (1. / (e * e)) * gg.transpose().kronecker(&gg);
    });
    mat.try_inverse()
        .unwrap()
        .diagonal()
        .iter()
        .map(|e2| (chi2 * e2).sqrt())
        .collect()
}
