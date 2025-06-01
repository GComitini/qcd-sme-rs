use lazy_static::lazy_static;
use log::info;
use peroxide::fuga::Algorithm;
use peroxide::util::plot::*;
use qcd_sme::R;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use t0_chempot::{alphastrong, gluon::propagator, init};
use t0_chempot::{Component, Config, BASEDIR, OMEPS};

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("t0fd_glueprop");
}

fn plot_propagators(config: &Config, corr: bool) {
    let corr_text = if corr { " (corrected)" } else { "" };
    let corr_ext = if corr { "corr_" } else { "" };

    let mus = config.fewchempots();
    let (mut data_t, mut data_l) = (Vec::with_capacity(mus.len()), Vec::with_capacity(mus.len()));

    for comp in &[Component::T, Component::L] {
        let comp_abbrev = comp.abbrev();
        let comp_lc = comp.to_string().to_lowercase();
        for mu in mus {
            let prop: Vec<R> = config
            .momenta()
            .par_iter()
            .map(|&p| {
                let prop = propagator(OMEPS, p, *mu, config, corr, *comp);
                info!("Computed{corr_text} {comp_abbrev} propagator at (mu, p) = ({mu:.4}, {p:.4}): D_{comp:?} = {prop:.5}");
                prop
            })
            .collect();

            match comp {
                Component::T => data_t.push(prop.clone()),
                Component::L => data_l.push(prop.clone()),
            };

            let mut outfile = BufWriter::new(
                std::fs::File::create(
                    THIS_BASEDIR.join(format!("{corr_ext}glueprop_{comp_lc}_mu_{mu:.4}.out")),
                )
                .unwrap_or_else(|e| panic!("Could not create output file: {e}")),
            );

            config.momenta().iter().zip(prop.iter()).for_each(|pp| {
                writeln!(outfile, "{}\t{}", pp.0, pp.1)
                    .unwrap_or_else(|e| panic!("Could not write to output file: {e}"));
            });

            let mut plot = Plot2D::new();
            plot.set_domain(config.momenta().clone());
            plot.insert_image(prop);
            plot.set_xlabel("$p$ [GeV]");
            plot.set_ylabel(&format!("$\\Delta_{{{comp:?}}}(p)$ [GeV$^{{-2}}$]"));
            plot.set_title(&format!("$\\mu = {mu:.4}$ GeV"));
            plot.set_ylim((0., config.ymax.unwrap_or(9.0)));
            plot.set_path(
                THIS_BASEDIR
                    .join(format!("{corr_ext}glueprop_{comp_lc}_mu_{mu:.4}.png"))
                    .to_str()
                    .unwrap(),
            );
            plot.savefig()
                .unwrap_or_else(|e| panic!("Could not save figure: {e}"));
        }
    }

    let legends: Vec<String> = mus.iter().map(|mu| format!("$\\mu={mu:.4}$ GeV")).collect();

    let mut plot = Plot2D::new();
    plot.set_domain(config.momenta().clone());
    data_t.into_iter().for_each(|d| {
        plot.insert_image(d);
    });
    plot.set_xlabel("$p$ [GeV]");
    plot.set_ylabel("$\\Delta_{{T}}(p)$ [GeV$^{{-2}}$]");
    plot.set_ylim((0., config.ymax.unwrap_or(9.0)));
    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.set_path(
        THIS_BASEDIR
            .join(format!("{corr_ext}glueprop_t.png"))
            .to_str()
            .unwrap(),
    );
    plot.savefig()
        .unwrap_or_else(|e| panic!("Could not save figure: {e}"));

    let mut plot = Plot2D::new();
    plot.set_domain(config.momenta().clone());
    data_l.into_iter().for_each(|d| {
        plot.insert_image(d);
    });
    plot.set_xlabel("$p$ [GeV]");
    plot.set_ylabel("$\\Delta_{{L}}(p)$ [GeV$^{{-2}}$]");
    plot.set_ylim((0., config.ymax.unwrap_or(9.0)));
    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.set_path(
        THIS_BASEDIR
            .join(format!("{corr_ext}glueprop_l.png"))
            .to_str()
            .unwrap(),
    );
    plot.savefig()
        .unwrap_or_else(|e| panic!("Could not save figure: {e}"));
}

fn main() {
    let config = init(Some(THIS_BASEDIR.as_path()));

    info!("Computing values of alpha_s...");
    let manychempots_init = config.manychempots_init();
    let (mumin, mumax) = (manychempots_init[0], manychempots_init[1]);
    let mut ass = Vec::with_capacity(8);
    for b in [true, false] {
        for c in [Component::T, Component::L] {
            for mm in [mumin, mumax] {
                ass.push(alphastrong(&config, mm, b, c))
            }
        }
    }
    let (asmin, asmax) = (ass.min(), ass.max());
    info!("alpha_s in [{asmin}, {asmax}]");

    plot_propagators(&config, false);
    if config.correctedfieldconfig().is_some() {
        plot_propagators(&config, true);
    }
}
