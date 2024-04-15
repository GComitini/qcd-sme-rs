use qcd_sme::ym::gluon::dressing_inv_landau;
use qcd_sme::ym::thermal::zero_momentum::gluon::dressing_t_inv_landau;
use qcd_sme::{C, I, R};

const TOL: R = 1E-12;
const MAX_ITER: usize = 1000000;
const Z0: C = C::new(-1., -1.);
const Z1: C = C::new(-0.5, -0.5);

fn find_root<F: Fn(C) -> C>(f: &F) -> Option<(C, C)> {
    qcd_sme::utils::find_root(f, (Z0, Z1), TOL, MAX_ITER)
}

fn find_pole<F: Fn(C) -> C>(header: &str, f: &F) {
    println!("{header}\n");
    let root = find_root(f);
    let root = match root {
        Some((root, value)) => {
            println!("Found root at s = {root} (f(s) = {value})");
            root
        }
        None => {
            println!("Found no root");
            return;
        }
    };
    let root = (-root).sqrt() * 0.656;
    println!(
        "For m = 0.656 GeV, this is equal to p = ({} + {} i) GeV.",
        root.re, root.im
    );
    println!();
}

fn find_thermal_pole<F: Fn(C) -> C>(header: &str, f: &F) {
    println!("{header}\n");
    let root = find_root(f);
    let root = match root {
        Some((root, value)) => {
            println!("Found root at omega = {root} (f(omega) = {value})");
            root
        }
        None => {
            println!("Found no root");
            return;
        }
    };
    let root = root * I;
    println!(
        "This is equal to epsilon = {} GeV, gamma = {} GeV.",
        root.re, -root.im
    );
    println!();
}

fn main() {
    find_pole("*** VACUUM POLE ***", &|z: C| {
        dressing_inv_landau(z, -0.876)
    });
    find_thermal_pole("*** T = 121 MeV POLE ***", &|z: C| {
        dressing_t_inv_landau(z, 0.656, 1. / 0.121, -0.836)
    });
    find_thermal_pole("*** T = 194 MeV POLE ***", &|z: C| {
        dressing_t_inv_landau(z, 0.550, 1. / 0.194, -0.696)
    });
    find_thermal_pole("*** T = 260 MeV POLE ***", &|z: C| {
        dressing_t_inv_landau(z, 0.450, 1. / 0.260, -0.416)
    });
    find_thermal_pole("*** T = 290 MeV POLE ***", &|z: C| {
        dressing_t_inv_landau(z, 0.450, 1. / 0.290, -0.476)
    });
    find_thermal_pole("*** T = 366 MeV POLE ***", &|z: C| {
        dressing_t_inv_landau(z, 0.450, 1. / 0.366, -0.196)
    });
    find_thermal_pole("*** T = 458 MeV POLE ***", &|z: C| {
        dressing_t_inv_landau(z, 0.450, 1. / 0.458, 0.214)
    });
}
