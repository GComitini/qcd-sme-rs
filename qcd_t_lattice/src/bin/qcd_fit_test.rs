use argmin::core::State;
use qcd_sme::consts::{set_default_max_iter_integral, set_default_tol_integral};
use qcd_t_lattice::fit::{param_errors, paramgradient3, ChiSquare3Params};
use qcd_t_lattice::gluon::propagator_l;
use qcd_t_lattice::{FieldConfig, OMEPS, R};

const F0: R = -0.876;
const MG: R = 0.656;

fn main() {
    set_default_max_iter_integral(30);
    set_default_tol_integral(1E-6);

    let fieldconfig = FieldConfig::new(3, MG, vec![(2, 0.35)]);

    let ndata = 40;
    let beta = 1. / 0.250;
    let data: Vec<(R, R, R)> = (0..=ndata)
        .map(|i| 0.01 + (i as R) / (ndata as R) * (4. - 0.01))
        .map(|p| (p, propagator_l(OMEPS, p, beta, F0, &fieldconfig).re, 0.5))
        .collect();

    let mq = 0.450;
    let ff = &|p, m, z, f0| {
        let fc = FieldConfig::new(3, m, vec![(2, mq)]);
        z * propagator_l(OMEPS, p, beta, f0, &fc).re
    };

    let cs3p = ChiSquare3Params { ff, data: &data };

    let exec = cs3p
        .executor(vec![0.3, 0.5, -0.5], true)
        .configure(|state| state.max_iters(100));
    let res = exec.run().unwrap();
    println!("{res}");

    let (chi2, params) = (
        res.state.get_best_cost(),
        res.state.get_best_param().unwrap(),
    );
    let (a, b, c) = (params[0], params[1], params[2]);
    let grad: Vec<Vec<R>> = data
        .iter()
        .map(|(p, _, _)| paramgradient3(ff, *p, a, b, c))
        .collect();

    let errors = param_errors(&data, &grad, chi2);
    println!("Errors: {errors:?}");
}
