use peroxide::fuga::{Plot, Plot2D};
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

fn find_pole<F: Fn(C) -> C>(header: &str, f: &F, res: &mut Vec<(C, R)>) {
    println!("{header}\n");
    let root = find_root(f);
    let root = match root {
        Some((root, value)) => {
            println!("Found root at s = {root} (f(s) = {value})");
            root
        }
        None => {
            panic!("Found no root");
        }
    };
    let residue = qcd_sme::utils::compute_residue(&|z| 1. / f(z), root) / root;
    let (z_r, z_i) = (residue.re, residue.im);
    let z_rat = z_i / z_r;
    let root = (-root).sqrt() * 0.656;
    res.push((root, z_rat));
    println!(
        "For m = 0.656 GeV, this is equal to p = ({} + {} i) GeV.",
        root.re, root.im
    );
    println!("The residue of the inverse is equal to {residue}, with im/re ratio = {z_rat}");
    println!();
}

fn find_thermal_pole<F: Fn(C) -> C>(header: &str, f: &F, res: &mut Vec<(C, R)>) {
    println!("{header}\n");
    let root = find_root(f);
    let root = match root {
        Some((root, value)) => {
            println!("Found root at omega = {root} (f(omega) = {value})");
            root
        }
        None => {
            panic!("Found no root");
        }
    };
    // 2*root is the conversion factor from a pole in p^2 to a pole in p
    // 1/root^2 gives us the residue of the propagator instead of that of the dressing function
    // -residue.im selects the pole on the same quadrant as that found by find_pole
    let residue = qcd_sme::utils::compute_residue(&|z| 1. / f(z), root) * 2. / root;
    let (z_r, z_i) = (residue.re, -residue.im);
    let z_rat = z_i / z_r;
    let root = root * I;
    res.push((root, z_rat));
    println!(
        "This is equal to epsilon = {} GeV, gamma = {} GeV.",
        root.re, -root.im
    );
    println!("The residue of the inverse is equal to {residue}, with im/re ratio = {z_rat}");
    println!();
}

fn main() {
    let mut res = Vec::new();
    let temperatures = [0., 0.121, 0.194, 0.260, 0.290, 0.366, 0.458];
    let masses = [0.656, 0.550, 0.450, 0.450, 0.450, 0.450];
    let f0s = [-0.836, -0.696, -0.416, -0.476, -0.196, 0.214];

    find_pole(
        "*** VACUUM POLE ***",
        &|z: C| dressing_inv_landau(z, -0.876),
        &mut res,
    );

    for (&t, (&m, &f0)) in temperatures
        .iter()
        .skip(1)
        .zip(masses.iter().zip(f0s.iter()))
    {
        find_thermal_pole(
            &format!("*** T = {} MeV POLE ***", (t * 1000.) as u32),
            &|z: C| dressing_t_inv_landau(z, m, 1. / t, f0),
            &mut res,
        );
    }

    let residues_ratio: Vec<R> = res.iter().map(|v: &(C, R)| v.1).collect();
    let residues_phase: Vec<R> = residues_ratio.iter().map(|rat| rat.atan()).collect();

    let mut plot = Plot2D::new();
    plot.set_domain(temperatures.to_vec());
    plot.insert_image(residues_ratio);
    plot.set_xlabel("$T$ (GeV)");
    plot.set_ylabel("Im$\\{R\\}$/Re$\\{R\\}$");
    plot.set_path("target/ym_gluon_poles_residue_ratio");
    plot.savefig().expect("Could not save figure");

    let mut plot = Plot2D::new();
    plot.set_domain(temperatures.to_vec());
    plot.insert_image(residues_phase);
    plot.set_xlabel("$T$ (GeV)");
    plot.set_ylabel("$\\vartheta(T)$");
    plot.set_path("target/ym_gluon_poles_residue_phase");
    plot.savefig().expect("Could not save figure");
}
