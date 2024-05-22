use peroxide::fuga::*;
use qcd_sme::consts::{get_number_of_colors, set_number_of_fermions};
use qcd_sme::qcd::gluon::dressing_crossed_inv_landau;
use std::{
    fs::File,
    io::{BufWriter, Result, Write},
};

static mut PLOT_INDEX: usize = 0;

fn do_dressing(
    momenta: &Vec<f64>,
    m: f64,
    mq: f64,
    f0: f64,
    z: f64,
    nf: u32,
    plot: &mut Plot2D,
) -> Result<()> {
    let mut outfile = BufWriter::new(File::create(&format!(
        "target/qcd_gluon_propagator_compare_{}.txt",
        unsafe { PLOT_INDEX }
    ))?);

    unsafe {
        PLOT_INDEX += 1;
    }

    set_number_of_fermions(nf);
    let nf_div_nc = (nf as f64) / (get_number_of_colors() as f64);

    let f0 = f0 + 2. * nf_div_nc / 27.;

    let m2 = m * m;

    let values = momenta
        .iter()
        .map(|p2| {
            let val = dressing_crossed_inv_landau(p2 / m2, mq / m, f0) / z;
            writeln!(outfile, "{}", format!("{p2}\t{val}")).ok();
            val
        })
        .collect();

    plot.insert_image(values);
    Ok(())
}

fn main() -> Result<()> {
    let (pmin, pmax) = (0.03, 10.);
    let dp = 0.01;
    let prange = pmax - pmin;
    let n = (prange / dp) as usize;

    let momenta: Vec<f64> = (0..=n)
        .map(|i| {
            let p = pmin + (i as f64) / (n as f64) * prange;
            p * p
        })
        .collect();

    let mut plot = Plot2D::new();
    plot.set_domain(momenta.clone());
    plot.set_path("target/qcd_gluon_propagator_compare.png");
    plot.set_xscale(PlotScale::Log);
    plot.set_yscale(PlotScale::Log);
    plot.set_xlim((0.001, 100.));
    plot.set_ylim((0.01, 100.));

    do_dressing(&momenta, 0.8, 0.65, -0.65, 0.4, 2, &mut plot)?;
    do_dressing(&momenta, 0.73, 0.65, -0.65, 0.4, 2, &mut plot)?;
    do_dressing(&momenta, 0.73, 0.65, -1.05, 1., 0, &mut plot)?;

    do_dressing(&momenta, 0.8, 0.5, -0.6, 0.4, 2, &mut plot)?;
    do_dressing(&momenta, 0.73, 0.5, -0.6, 0.4, 2, &mut plot)?;
    do_dressing(&momenta, 0.73, 0.5, -1.05, 1., 0, &mut plot)?;

    plot.savefig().ok();
    Ok(())
}
