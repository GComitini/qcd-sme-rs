use qcd_sme::consts::{set_default_max_iter_integral, set_default_tol_integral};
use qcd_sme::ym::gluon::propagator_landau;
use qcd_t_lattice::fit::{ChiSquare2Params, ChiSquare3Params, ChiSquare4Params};
use qcd_t_lattice::{OMEPS, R};

fn main() {
    set_default_max_iter_integral(30);
    set_default_tol_integral(1E-6);

    let ndata = 40;
    let data: Vec<(R, R, R)> = (0..=ndata)
        .map(|i| 0.01 + (i as R) / (ndata as R) * (4. - 0.01))
        .map(|p| {
            let m = 0.656;
            let s = p / m;
            let s = s * s;
            (p, 2.6 * propagator_landau(s, -0.876) / (m * m), 0.5)
        })
        .collect();

    let cs2p = ChiSquare2Params {
        ff: &|p, m, z| {
            let s = p / m;
            let s = s * s;
            z * propagator_landau(s, -0.876) / (m * m)
        },
        data: &data,
    };

    let exec = cs2p.executor(vec![0.5, 1.], true);
    let res = exec.run().unwrap();
    println!("{res}");

    let cs3p = ChiSquare3Params {
        ff: &|p, m, z, f0| {
            let s = p / m;
            let s = s * s;
            z * propagator_landau(s, f0) / (m * m)
        },
        data: &data,
    };

    let exec = cs3p.executor(vec![0.5, 1., -0.5], true);
    let res = exec.run().unwrap();
    println!("{res}");

    let cs4p = ChiSquare4Params {
        ff: &|p, m, z, f0, y| {
            let s = (p + y) / m;
            let s = s * s;
            z * propagator_landau(s, f0) / (m * m)
        },
        data: &data,
    };

    let exec = cs4p.executor(vec![0.5, 1., -0.5, 0.3], true);
    let res = exec.run().unwrap();
    println!("{res}");

    let cs2p_ft = ChiSquare2Params {
        ff: &|p, m, z| {
            z * qcd_sme::ym::thermal::gluon::propagator_t_landau(OMEPS, p, m, 1. / 0.01, -0.876).re
        },
        data: &data,
    };

    let exec = cs2p_ft
        .executor(vec![0.5, 1.], true)
        .configure(|state| state.target_cost(1E-6));
    let res = exec.run().unwrap();
    println!("{res}");

    let cs3p_ft = ChiSquare3Params {
        ff: &|p, m, z, f0| {
            z * qcd_sme::ym::thermal::gluon::propagator_t_landau(OMEPS, p, m, 1. / 0.01, f0).re
        },
        data: &data,
    };

    let exec = cs3p_ft
        .executor(vec![0.5, 1., -0.5], true)
        .configure(|state| state.target_cost(1E-6));
    let res = exec.run().unwrap();
    println!("{res}");

    let cs3p_ft = ChiSquare3Params {
        ff: &|p, m, z, f0| {
            z * qcd_sme::ym::thermal::gluon::propagator_t_landau(OMEPS, p, m, 1. / 0.15, f0).re
        },
        data: &data,
    };

    let exec = cs3p_ft
        .executor(vec![0.5, 1., -0.5], true)
        .configure(|state| state.target_cost(1E-6));
    let res = exec.run().unwrap();
    println!("{res}");
}
