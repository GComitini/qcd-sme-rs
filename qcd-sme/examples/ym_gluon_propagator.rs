// Peroxide is used to plot the propagators (and internally for the thermal
// integrals)
use peroxide::fuga::*;
use qcd_sme::ym::thermal::gluon::{propagator_l_landau, propagator_t_landau};

fn main() {
    /* DEFINITIONS: MOMENTA, TEMPERATURE, ETC. */
    let (pmin, pmax) = (0.2, 3.);
    let prange = pmax - pmin;
    let dp = 0.01;
    let n = (prange / dp) as usize;

    let momenta: Vec<f64> = (0..=n)
        .map(|i| pmin + prange * (i as f64) / (n as f64))
        .collect();

    let t = 0.260;
    let beta = 1. / t;
    let renpoint = 2.;
    let om = 0.001;

    /* 3DIMENSIONALLY-TRANSVERSE PROPAGATOR  */
    println!("*** TRANSVERSE PROJECTION ***");
    let m = 0.450;
    let f0 = -0.42;

    let norm = propagator_t_landau(om, renpoint, m, beta, f0).re * (renpoint * renpoint);

    let glplt = momenta
        .iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_t_landau(om, p, m, beta, f0).re / norm;
            eprintln!("Computed {res}.");
            res
        })
        .collect();

    /* 3DIMENSIONALLY-LONGITUDINAL PROPAGATOR  */
    println!("\n*** LONGITUDINAL PROJECTION ***");
    let m = 0.425;
    let f0 = -1.42;

    let norm = propagator_l_landau(om, renpoint, m, beta, f0).re * (renpoint * renpoint);

    let glpll = momenta
        .iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_l_landau(om, p, m, beta, f0).re / norm;
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
