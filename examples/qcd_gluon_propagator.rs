// Peroxide is used to plot the propagators (and internally for the thermal
// integrals)
use peroxide::fuga::*;
use qcd_sme::consts::get_matsubara_reg;
use qcd_sme::qcd::thermal::gluon::{
    propagator_l_landau, propagator_l_zero_temp_landau, propagator_t_landau,
    propagator_t_zero_temp_landau,
};
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
    let mu = 1.2;
    let m = 0.656;
    let f0 = -0.876;
    let om = get_matsubara_reg();

    let temperatures = [0.100, 0.075, 0.050, 0.025, 0.010, 0.005];

    /* 3DIMENSIONALLY-TRANSVERSE POLARIZATION  */
    let mut plot = peroxide::util::plot::Plot2D::new();
    plot.set_domain(momenta.clone());
    plot.set_path(&format!("target/qcd_gluon_propagator_transverse.png"));
    let mut legends = Vec::new();

    for &t in &temperatures {
        let beta = 1. / t;
        let glplt = momenta
            .par_iter()
            .map(|&p| {
                eprint!("Computing {p}... ");
                let res = propagator_t_landau(om, p, m, beta, mu, f0).re;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(glplt);
        legends.push(format!("T = {t} GeV",));
    }

    let glplt = momenta
        .par_iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_t_zero_temp_landau(om, p, m, mu, f0).re;
            eprintln!("Computed {res}.");
            res
        })
        .collect();

    plot.insert_image(glplt);
    legends.push(format!("T = 0 GeV",));

    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.savefig().expect("Could not save figure");

    /* 3DIMENSIONALLY-LONGITUDINAL POLARIZATION  */
    let mut plot = peroxide::util::plot::Plot2D::new();
    plot.set_domain(momenta.clone());
    plot.set_path(&format!("target/qcd_gluon_propagator_longitudinal.png"));
    let mut legends = Vec::new();

    for &t in &temperatures {
        let beta = 1. / t;
        let glplt = momenta
            .par_iter()
            .map(|&p| {
                eprint!("Computing {p}... ");
                let res = propagator_l_landau(om, p, m, beta, mu, f0).re;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(glplt);
        legends.push(format!("T = {t} GeV",));
    }

    let glplt = momenta
        .par_iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_l_zero_temp_landau(om, p, m, mu, f0).re;
            eprintln!("Computed {res}.");
            res
        })
        .collect();

    plot.insert_image(glplt);
    legends.push(format!("T = 0 GeV",));

    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.savefig().expect("Could not save figure");
}
