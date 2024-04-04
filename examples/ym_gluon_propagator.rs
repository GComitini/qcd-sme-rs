// Peroxide is used to plot the propagators (and internally for the thermal
// integrals)
use peroxide::fuga::*;
use qcd_sme::ym::thermal::gluon::{propagator_l_landau, propagator_t_landau};
// Rayon is used to parallelize the computation of the values of the propagator
// using par_iter
use rayon::prelude::*;

fn main() {
    /* DEFINITIONS: MOMENTA, TEMPERATURE, ETC. */
    let (pmin, pmax) = (0.1, 4.);
    let prange = pmax - pmin;
    let dp = 0.01;
    let n = (prange / dp) as usize;

    let momenta: Vec<f64> = (0..=n)
        .map(|i| pmin + prange * (i as f64) / (n as f64))
        .collect();

    let t = 0.260;
    let beta = 1. / t;
    let renpoint = 4.;

    /* 3DIMENSIONALLY-TRANSVERSE PROPAGATOR  */
    let m = 0.450;
    let f0 = -0.42;

    let norm = propagator_t_landau(0.001, renpoint, m, beta, f0).re * (renpoint * renpoint);

    let glplt = momenta
        .par_iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_t_landau(0.001, p, m, beta, f0).re / norm;
            eprintln!("Computed {res}.");
            res
        })
        .collect();

    /* 3DIMENSIONALLY-LONGITUDINAL PROPAGATOR  */
    let m = 0.425;
    let f0 = -1.42;

    let norm = propagator_l_landau(0.001, renpoint, m, beta, f0).re * (renpoint * renpoint);

    let glpll = momenta
        .par_iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_l_landau(0.001, p, m, beta, f0).re / norm;
            eprintln!("Computed {res}.");
            res
        })
        .collect();

    /* PLOT */
    let mut plot = peroxide::util::plot::Plot2D::new();
    plot.set_domain(momenta);
    plot.insert_image(glplt);
    plot.insert_image(glpll);
    plot.set_path("target/ym_gluon_propagator.png");
    plot.savefig().expect("Could not save figure");
}
