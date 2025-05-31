use qcd_sme::consts::{set_default_max_iter_integral, set_default_tol_integral};
use qcd_t_lattice::{
    fit::{readfile, ChiSquare3Params},
    OMEPS, R,
};

fn main() {
    set_default_max_iter_integral(25);
    set_default_tol_integral(1E-5);

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
        .expect("You need to pass the temperature as argument!")
        .parse::<R>()
        .expect("The temperature argument must be a float!");

    let cs3p_ft = ChiSquare3Params {
        ff: &|p, m, z, f0| match comp.as_str() {
            "l" | "L" => {
                z * qcd_sme::ym::thermal::gluon::propagator_l_landau(OMEPS, p, m, 1. / t, f0).re
            }
            "t" | "T" => {
                z * qcd_sme::ym::thermal::gluon::propagator_t_landau(OMEPS, p, m, 1. / t, f0).re
            }
            _ => unreachable!(),
        },
        data: &data,
    };

    let exec = cs3p_ft
        .executor(vec![0.65, 2., -0.8], true)
        .configure(|state| state.target_cost(1E-6).max_iters(100));
    let res = exec.run().unwrap();
    println!("{res}");
}
