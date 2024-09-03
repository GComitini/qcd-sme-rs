use peroxide::fuga::{Markers, Plot, Plot2D, PlotType};
use qcd_sme::ym::gluon::dressing_inv_landau;
use qcd_sme::ym::thermal::zero_momentum::gluon::dressing_t_inv_landau;
use qcd_sme::{Num, C, I, R};
use std::io::{BufWriter, Write};

const MG: R = 0.656;

const TOL: R = 1E-12;
const MAX_ITER: usize = 1000000;
const Z0: C = C::new(-1., -1.);
const Z1: C = C::new(-0.5, -0.5);

fn find_root<F: Fn(C) -> C>(f: &F) -> Option<(C, C)> {
    qcd_sme::utils::find_root(f, (Z0, Z1), TOL, MAX_ITER)
}

// This actually finds zeroes
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
    let root = (-root).sqrt() * MG;
    res.push((root, z_rat));
    println!(
        "For m = {MG:.3} GeV, this is equal to p = ({} + {} i) GeV.",
        root.re, root.im
    );
    println!("The residue of the inverse is equal to {residue}, with im/re ratio = {z_rat}");
    println!();
}

// This actually finds zeroes
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

    // Check that there are four symmetric roots
    if (f(-root)).abs() > TOL || (f(root.conj())).abs() > TOL || (f(-root.conj())).abs() > TOL {
        panic!("Roots are not symmetric!");
    }

    // +/- 2*root is the conversion factor from a pole in p to a pole in p^2
    // 1/root^2 gives us the residue of the propagator instead of that of the dressing function
    // -residue.im and conj() select the pole on the first quadrant
    let residue = qcd_sme::utils::compute_residue(&|z| 1. / f(z), root) * 2. / root;
    let residue = residue.conj();
    let (z_r, z_i) = (residue.re, residue.im);
    let z_rat = z_i / z_r;
    let root = (root * I).conj();
    res.push((root, z_rat));
    println!(
        "This is equal to epsilon = {} GeV, gamma = {} GeV.",
        root.re, root.im
    );
    println!("The residue of the inverse is equal to {residue}, with im/re ratio = {z_rat}");
    println!();
}

fn main() {
    let mut res = Vec::new();

    let f0 = -0.876;

    let temperatures = [0., 0.121, 0.194, 0.260, 0.290, 0.366, 0.458];
    let masses = [0.675, 0.725, 0.775, 0.725, 0.800, 0.900];
    let df0s = [0.05, 0.10, 0.30, 0.40, 0.50, 0.60];

    find_pole(
        "*** VACUUM POLE ***",
        &|z: C| dressing_inv_landau(z, f0),
        &mut res,
    );

    for (&t, (&m, &df0)) in temperatures
        .iter()
        .skip(1)
        .zip(masses.iter().zip(df0s.iter()))
    {
        find_thermal_pole(
            &format!("*** T = {} MeV POLE ***", (t * 1000.) as u32),
            &|z: C| dressing_t_inv_landau(z, m, 1. / t, df0 + f0),
            &mut res,
        );
    }

    let poles_re: Vec<R> = res.iter().map(|v: &(C, R)| v.0.re).collect();
    let poles_im: Vec<R> = res.iter().map(|v: &(C, R)| v.0.im).collect();
    let residues_ratio: Vec<R> = res.iter().map(|v: &(C, R)| v.1).collect();
    let residues_phase: Vec<R> = residues_ratio.iter().map(|rat| rat.atan()).collect();

    let out_dir = std::path::Path::new("target/ym_thermal_gluon_poles");
    std::fs::create_dir_all(out_dir)
        .unwrap_or_else(|_| panic!("Could not create {}", out_dir.to_string_lossy()));

    let out_filename = out_dir.join("poles.out");
    let mut out_file = BufWriter::new(
        std::fs::File::create(&out_filename)
            .unwrap_or_else(|_| panic!("Could not create {}", out_filename.to_string_lossy())),
    );

    temperatures.iter().zip(res.iter()).for_each(|(t, val)| {
        writeln!(out_file, "{t:.3}\t{}\t{}", val.0, val.1)
            .unwrap_or_else(|_| panic!("Could not write to {}", out_filename.to_string_lossy()))
    });

    let mut plot = Plot2D::new();
    plot.set_domain(temperatures.to_vec());
    plot.set_plot_type(vec![(0, PlotType::Scatter), (1, PlotType::Scatter)]);
    plot.set_marker(vec![(0, Markers::Square), (1, Markers::Square)]);
    plot.insert_image(poles_re);
    plot.insert_image(poles_im);
    plot.set_xlabel("$T$ (GeV)");
    plot.set_ylabel("Re$\\{p\\}$, Im$\\{p\\}$ (GeV)");
    plot.set_legend(vec!["Re$\\{p\\}$", "Im$\\{p\\}$"]);
    plot.set_ylim((0., 1.));
    plot.set_path(&out_dir.join("poles").to_string_lossy());
    plot.savefig().expect("Could not save figure");

    let mut plot = Plot2D::new();
    plot.set_domain(temperatures.to_vec());
    plot.set_plot_type(vec![(0, PlotType::Scatter)]);
    plot.set_marker(vec![(0, Markers::Square)]);
    plot.insert_image(residues_ratio);
    plot.set_xlabel("$T$ (GeV)");
    plot.set_ylabel("Im$\\{R\\}$/Re$\\{R\\}$");
    plot.set_path(&out_dir.join("residue_ratio").to_string_lossy());
    plot.savefig().expect("Could not save figure");

    let mut plot = Plot2D::new();
    plot.set_domain(temperatures.to_vec());
    plot.set_plot_type(vec![(0, PlotType::Scatter)]);
    plot.set_marker(vec![(0, Markers::Square)]);
    plot.insert_image(residues_phase);
    plot.set_xlabel("$T$ (GeV)");
    plot.set_ylabel("$\\vartheta(T)$");
    plot.set_path(&out_dir.join("residue_phase").to_string_lossy());
    plot.savefig().expect("Could not save figure");
}
