use peroxide::{fuga::*, util::plot};
use qcd_sme::consts::set_default_tol_integral;
use qcd_sme::qcd::thermal::quark::{chemical_potential, chemical_potential_zero_temp};
use rayon::prelude::*;

fn main() {
    let (nmin, nmax) = (-2., 2.);
    let nrange = nmax - nmin;
    let dn = 0.001;
    let k = (nrange / dn) as u32;
    let densities: Vec<f64> = (0..k)
        .map(|i| nmin + nrange * (i as f64) / (k as f64))
        .collect();
    let mq = 0.5;
    let c: usize = 1250;
    let detail_range = c..(k as usize - c);
    set_default_tol_integral(1E-5);

    let mut plot = plot::Plot2D::new();
    plot.set_domain(densities.clone());
    plot.set_path(&"target/qcd_quark_chemical_potential.png".to_string());
    let mut plot_detail = plot::Plot2D::new();
    plot_detail.set_domain(densities[detail_range.clone()].to_vec());
    plot_detail.set_path(&"target/qcd_quark_chemical_potential_detail.png".to_string());
    let mut legends = Vec::new();

    let temperatures = [
        0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.,
    ];

    legends.push("T = 0.0 GeV".to_string());
    let vals: Vec<f64> = densities
        .par_iter()
        .map(|&n| {
            eprintln!("Computing T = 0.0 GeV, mu = {n} GeV");
            chemical_potential_zero_temp(mq, n)
        })
        .collect();
    plot_detail.insert_image(vals[detail_range.clone()].to_vec());
    plot.insert_image(vals);

    for t in temperatures {
        legends.push(format!("T = {t} GeV"));
        let vals: Vec<f64> = densities
            .par_iter()
            .map(|&mu| {
                eprintln!("Computing T = {t} GeV, mu = {mu} GeV");
                chemical_potential(mq, 1. / t, mu)
            })
            .collect();
        if t <= 0.2 {
            plot_detail.insert_image(vals[detail_range.clone()].to_vec());
        }
        plot.insert_image(vals);
    }

    plot.set_legend(legends.iter().map(|l| l.as_str()).collect());
    plot.savefig().expect("Could not save figure");

    plot_detail.set_legend(legends.iter().map(|l| l.as_str()).collect());
    plot_detail.savefig().expect("Could not save figure");
}
