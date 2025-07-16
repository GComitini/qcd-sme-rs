use peroxide::prelude::*;
use qcd_sme::low_level::oneloop::thermal::gluon::zero_momentum::polarization_quark_l_thermal_part_landau_i as polarization_quark_l_thermal_part_landau_i_zero_momentum;
use qcd_sme::low_level::oneloop::thermal::gluon::{
    polarization_quark_l_thermal_part_landau_i, polarization_quark_t_thermal_part_landau_i,
};
use rayon::prelude::*;

fn main() {
    /* DEFINITIONS: MOMENTA, TEMPERATURE, ETC. */
    let (qmin, qmax) = (0.01, 1.);
    let qrange = qmax - qmin;
    let dq = 0.001;
    let n = (qrange / dq) as usize;

    let momenta: Vec<f64> = (0..=n)
        .map(|i| qmin + qrange * (i as f64) / (n as f64))
        .collect();

    /* PROPAGATOR PARAMETERS */
    let t = 0.260;
    let beta = 1. / t;
    let mq = 0.300;
    let om = 0.1;

    /* PLOT INIT */
    let mut plot = peroxide::util::plot::Plot2D::new();
    plot.set_domain(momenta.clone());
    plot.set_path("target/gluon_quark_polarization.png");
    plot.set_xlabel("q (GeV) - integration variable");

    let mut legend = Vec::new();

    /* CALCULATIONS */

    for p in [0.1, 0.05, 0.01] {
        let glplt = momenta
            .par_iter()
            .map(|&q| {
                eprint!("Computing {q} (longitudinal)... ");
                let res = polarization_quark_l_thermal_part_landau_i(q, om, p, mq, beta, 1.2).re;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(glplt);
        legend.push(format!("p = {p} GeV (long.)"));
        let glplt = momenta
            .par_iter()
            .map(|&q| {
                eprint!("Computing {q} (transverse)... ");
                let res = polarization_quark_t_thermal_part_landau_i(q, om, p, mq, beta, 1.2).re;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(glplt);
        legend.push(format!("p = {p} GeV (trans.)"));
    }

    let glplt = momenta
        .par_iter()
        .map(|&q| {
            eprint!("Computing {q}... ");
            let res =
                polarization_quark_l_thermal_part_landau_i_zero_momentum(q, 0.1, mq, beta, 1.2);
            eprintln!("Computed {res}.");
            res
        })
        .collect();
    plot.insert_image(glplt);
    legend.push("p = 0 GeV".to_string());

    /* PLOT */

    plot.set_legend(legend.iter().map(|s| s.as_str()).collect());
    plot.savefig().expect("WTF");
}
