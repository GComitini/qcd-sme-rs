use peroxide::{fuga::*, util::plot};
use qcd_sme::qcd::thermal::quark::{charge_density, charge_density_zero_temp};
use rayon::prelude::*;

fn main() {
    let (mumin, mumax) = (-2., 2.);
    let murange = mumax - mumin;
    let dmu = 0.001;
    let n = (murange / dmu) as u32;
    let chempots: Vec<f64> = (0..n)
        .map(|i| mumin + murange * (i as f64) / (n as f64))
        .collect();
    let mq = 0.5;
    let c: usize = 1250;
    let detail_range = c..(n as usize - c);

    let mut plot = plot::Plot2D::new();
    plot.set_domain(chempots.clone());
    plot.set_path(&"target/qcd_quark_charge_density.png".to_string());
    let mut plot_detail = plot::Plot2D::new();
    plot_detail.set_domain(chempots[detail_range.clone()].to_vec());
    plot_detail.set_path(&"target/qcd_quark_charge_density_detail.png".to_string());
    let mut legends = Vec::new();

    let temperatures = [
        0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.,
    ];

    legends.push("T = 0.0 GeV".to_string());
    let vals: Vec<f64> = chempots
        .par_iter()
        .map(|&mu| {
            eprintln!("Computing T = 0.0 GeV, mu = {mu} GeV");
            charge_density_zero_temp(mq, mu)
        })
        .collect();
    plot_detail.insert_image(vals[detail_range.clone()].to_vec());
    plot.insert_image(vals);

    for t in temperatures {
        legends.push(format!("T = {t} GeV"));
        let vals: Vec<f64> = chempots
            .par_iter()
            .map(|&mu| {
                eprintln!("Computing T = {t} GeV, mu = {mu} GeV");
                charge_density(mq, 1. / t, mu)
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
