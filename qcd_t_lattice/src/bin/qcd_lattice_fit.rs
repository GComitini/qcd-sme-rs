use argmin::core::State;
use qcd_sme::consts::{set_default_max_iter_integral, set_default_tol_integral};
use qcd_t_lattice::{
    fit::{
        param_errors, paramgradient2, paramgradient3, readfile, ChiSquare3Params, FitType,
        RenChiSquare2Params,
    },
    gluon::*,
    FieldConfig, OMEPS, R,
};

const P0: R = 1E-2;
const PREN: R = 4.;

fn main() {
    let data = readfile("data.txt");

    let comp = std::env::args()
        .nth(1)
        .expect("You need to pass the component as argument!");

    match comp.as_str() {
        "l" | "L" | "t" | "T" => {}
        _ => panic!("The component argument can only be t, T, l or L!"),
    }

    let t = std::env::args()
        .nth(2)
        .expect("You need to pass the temperature in GeV as argument!")
        .parse::<R>()
        .expect("The temperature argument must be a float!");

    let mq = std::env::args()
        .nth(3)
        .expect("You need to pass the quark mass in GeV as argument!")
        .parse::<R>()
        .expect("The quark mass argument must be a float!");

    let niter = std::env::args()
        .nth(4)
        .expect("You need to pass the number of iterations as argument!")
        .parse::<u64>()
        .expect("The number of iterations argument must be an integer!");

    let (mdef, zdef, f0def) = (0.7, 2., -0.8);

    let (fittype, m, z, f0) = match std::env::args()
        .nth(5)
        .expect("You need to pass the type of fit as argument!")
        .as_str()
    {
        "ren" => {
            if let (Some(m), Some(f0)) = (std::env::args().nth(6), std::env::args().nth(7)) {
                (
                    FitType::RenormZ,
                    m.replace(",", "")
                        .parse::<R>()
                        .expect("The gluon mass argument must be a float!"),
                    zdef,
                    f0.parse::<R>().expect("The f0 argument must be a float!"),
                )
            } else {
                println!("WARNING: using default parameters");
                (FitType::RenormZ, mdef, zdef, f0def)
            }
        }
        "unren" => {
            if let (Some(m), Some(z), Some(f0)) = (
                std::env::args().nth(6),
                std::env::args().nth(7),
                std::env::args().nth(8),
            ) {
                (
                    FitType::FittedZ,
                    m.replace(",", "")
                        .parse::<R>()
                        .expect("The gluon mass argument must be a float!"),
                    z.replace(",", "")
                        .parse::<R>()
                        .expect("The z argument must be a float!"),
                    f0.parse::<R>().expect("The f0 argument must be a float!"),
                )
            } else {
                println!("WARNING: using default parameters");
                (FitType::FittedZ, mdef, zdef, f0def)
            }
        }
        el => panic!("Unknown renormalization type {el}"),
    };

    let nextarg = if fittype == FitType::RenormZ { 8 } else { 9 };

    let mii = std::env::args().nth(nextarg).map_or(30, |mii| {
        mii.parse::<u32>()
            .expect("The max integral iterations argument must be an integer!")
    });

    set_default_max_iter_integral(mii);
    set_default_tol_integral(1E-7);

    match fittype {
        FitType::FittedZ => {
            let stp = vec![m, z, f0];

            let ff = &|p, m, z, f0| {
                let p = if p == 0. { P0 } else { p };
                let fieldconfig = FieldConfig::new(3, m, vec![(2, mq)]);
                match comp.as_str() {
                    "l" | "L" => z * propagator_l(OMEPS, p, 1. / t, f0, &fieldconfig).re,
                    "t" | "T" => z * propagator_t(OMEPS, p, 1. / t, f0, &fieldconfig).re,
                    _ => unreachable!(),
                }
            };

            let cs3p_ft = ChiSquare3Params { ff, data: &data };

            let exec = cs3p_ft
                .executor(stp.clone(), true)
                .configure(|state| state.max_iters(niter));
            let res = exec.run().unwrap();

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

            println!("----------");
            println!(
        "Component: {}. T: {t} GeV. mq: {mq}. niter: {niter}. Starting point: {stp:?}. Max integral iterations: {mii}.\n",
        comp.to_uppercase()
    );
            println!("{res}");

            println!("Errors: {errors:?}");
            println!("----------");
        }
        FitType::RenormZ => {
            let stp = vec![m, f0];

            let ff = &|p, m, f0| {
                let p = if p == 0. { P0 } else { p };
                let fieldconfig = FieldConfig::new(3, m, vec![(2, mq)]);
                match comp.as_str() {
                    "l" | "L" => propagator_l(OMEPS, p, 1. / t, f0, &fieldconfig).re,
                    "t" | "T" => propagator_t(OMEPS, p, 1. / t, f0, &fieldconfig).re,
                    _ => unreachable!(),
                }
            };

            let cs2p_ft = RenChiSquare2Params {
                ff,
                data: &data,
                pren: PREN,
            };

            let exec = cs2p_ft
                .executor(stp.clone(), true)
                .configure(|state| state.max_iters(niter));
            let res = exec.run().unwrap();

            let (chi2, params) = (
                res.state.get_best_cost(),
                res.state.get_best_param().unwrap(),
            );
            let (a, b) = (params[0], params[1]);
            let grad: Vec<Vec<R>> = data
                .iter()
                .map(|(p, _, _)| paramgradient2(ff, *p, a, b))
                .collect();

            let errors = param_errors(&data, &grad, chi2);

            println!("----------");
            println!(
        "Component: {}. T: {t} GeV. mq: {mq}. niter: {niter}. Starting point: {stp:?}. Max integral iterations: {mii}.\n",
        comp.to_uppercase()
    );
            println!("{res}");

            println!("Errors: {errors:?}");
            println!("----------");
        }
    };
}
