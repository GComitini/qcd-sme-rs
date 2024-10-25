use lazy_static::lazy_static;
use log::info;
use peroxide::prelude::*;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use t0_chempot::{gluon::propagator, init};
use t0_chempot::{Component, Config, BASEDIR, OMEPS, PBASE_DEFAULT, R};

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("t0fd_gluemass");
}

// This is defined as the square root of the inverse propagator at p = 0
fn mass(mu: R, config: &Config, corr: bool, comp: Component) -> R {
    let m = (1. / propagator(OMEPS, PBASE_DEFAULT, mu, config, corr, comp)).sqrt();
    info!(
        "Computed{} T mass for mu = {mu:.4}: m_T = {m:.4}",
        if corr { " (corrected)" } else { "" }
    );
    m
}

fn plot_masses(config: &Config, corr: bool) {
    let corr_text = if corr { " (corrected)" } else { "" };
    let corr_ext = if corr { "corr_" } else { "" };

    let mass_t: Vec<R> = config
        .manychempots()
        .par_iter()
        .map(|mu| mass(*mu, config, corr, Component::T))
        .collect();

    let mass_l: Vec<R> = config
        .manychempots()
        .par_iter()
        .map(|mu| mass(*mu, config, corr, Component::L))
        .collect();

    let mut outfile = BufWriter::new(
        std::fs::File::create(THIS_BASEDIR.join(format!(
            "{corr_ext}gluemass_f0_{:.4}_pren_{:.4}.out",
            config.f0_init(),
            config.pren
        )))
        .unwrap_or_else(|e| panic!("Could not create output file: {e}")),
    );

    config
        .manychempots()
        .iter()
        .zip(mass_t.iter().zip(mass_l.iter()))
        .for_each(|(mu, (mt, ml))| {
            writeln!(outfile, "{mu:.4}\t{mt:.4}\t{ml:.4}")
                .unwrap_or_else(|e| panic!("Could not write to output file: {e}"));
        });

    let mut plot = Plot2D::new();
    plot.set_domain(config.manychempots().clone());
    plot.insert_image(mass_t);
    plot.insert_image(mass_l);
    plot.set_ylim((0., 1.5));
    plot.set_title(&format!("$T=0$ gluon (propagator) mass{corr_text}"));
    plot.set_path(
        THIS_BASEDIR
            .join(format!(
                "{corr_ext}gluemass_f0_{:.4}_pren_{:.4}.png",
                config.f0_init(),
                config.pren
            ))
            .to_str()
            .unwrap(),
    );
    plot.set_legend(vec!["Transverse", "Longitudinal"]);
    plot.savefig()
        .unwrap_or_else(|e| panic!("Could not save figure: {e}"));
}

fn main() {
    let config = init(Some(THIS_BASEDIR.as_path()));
    plot_masses(&config, false);
    if config.correctedfieldconfig().is_some() {
        plot_masses(&config, true);
    }
}
