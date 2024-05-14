// Peroxide is used to plot the propagators (and internally for the thermal
// integrals)
use peroxide::fuga::*;
use qcd_sme::ym::gluon::propagator_landau;
use qcd_sme::ym::thermal::gluon::{propagator_l_landau, propagator_t_landau};
use qcd_sme::R;
use std::fs;

fn main() {
    /* DEFINITIONS: MOMENTA, TEMPERATURE, ETC. */
    let (pmin, pmax) = (0.2, 3.);
    let prange = pmax - pmin;
    let dp = 0.01;
    let n = (prange / dp) as usize;
    let momenta: Vec<R> = (0..=n).map(|i| pmin + dp * (i as R)).collect();
    let renpoint = 4.;
    let om = 0.001;

    let temperatures = [0.121, 0.194, 0.260, 0.290, 0.366, 0.458];
    let ms_l = [0.550, 0.425, 0.425, 0.275, 0.150];
    let f0s_l = [-0.886, -1.099, -1.421, -0.966, -0.596];
    let ms_t = [0.656, 0.550, 0.450, 0.450, 0.450, 0.450];
    let f0s_t = [-0.836, -0.696, -0.416, -0.476, -0.196, 0.214];

    let targetdir = std::path::Path::new("target/ym_gluon_propagator");
    fs::create_dir_all(targetdir).expect(&format!(
        "Could not crate directory {}",
        targetdir.to_string_lossy()
    ));

    /* VACUUM PROPAGATOR */
    let norm =
        propagator_landau((renpoint * renpoint) / (0.656 * 0.656), -0.876) * (renpoint * renpoint);
    let vacuum_vals: Vec<R> = momenta
        .iter()
        .map(|&p| {
            eprint!("Computing {p}... ");
            let res = propagator_landau((p * p) / (0.656 * 0.656), -0.876) / norm;
            eprintln!("Computed {res}.");
            res
        })
        .collect();

    /* SUBCRITICAL TEMPERATURE, LONGITUDINAL PROPAGATOR */
    let mut plot = Plot2D::new();
    plot.set_xlabel("$p$ [GeV]");
    plot.set_ylabel("$\\Delta_{L}(p)$ [GeV$^{-2}$]");
    plot.set_path(
        &targetdir
            .join("gluon_propagator_l_low_t.png")
            .to_string_lossy(),
    );
    plot.set_domain(momenta.clone());

    plot.insert_image(vacuum_vals.clone());
    let mut legends = vec![String::from("$T=0$ MeV")];

    for i in 0..=2 {
        let t = temperatures[i];
        let beta = 1. / t;
        let m = ms_l[i];
        let f0 = f0s_l[i];
        let norm = propagator_l_landau(om, renpoint, m, beta, f0).re * (renpoint * renpoint);
        let vals = momenta
            .iter()
            .map(|&p| {
                eprint!("Computing {p}... ");
                let res = propagator_l_landau(om, p, m, beta, f0).re / norm;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(vals);
        legends.push(format!("$T={}$ MeV", (t * 1000.) as u32));
    }
    plot.set_legend(legends.iter().map(|l| l.as_str()).collect());
    plot.savefig().expect("Could not save figure");

    /* SUPERCRITICAL TEMPERATURE, LONGITUDINAL PROPAGATOR */
    let mut plot = Plot2D::new();
    plot.set_xlabel("$p$ [GeV]");
    plot.set_ylabel("$\\Delta_{L}(p)$ [GeV$^{-2}$]");
    plot.set_path(
        &targetdir
            .join("gluon_propagator_l_high_t.png")
            .to_string_lossy(),
    );
    plot.set_domain(momenta.clone());
    let mut legends = vec![];
    for i in 3..=4 {
        let t = temperatures[i];
        let beta = 1. / t;
        let m = ms_l[i];
        let f0 = f0s_l[i];
        let norm = propagator_l_landau(om, renpoint, m, beta, f0).re * (renpoint * renpoint);
        let vals = momenta
            .iter()
            .map(|&p| {
                eprint!("Computing {p}... ");
                let res = propagator_l_landau(om, p, m, beta, f0).re / norm;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(vals);
        legends.push(format!("$T={}$ MeV", (t * 1000.) as u32));
    }
    plot.set_legend(legends.iter().map(|l| l.as_str()).collect());
    plot.savefig().expect("Could not save figure");

    /* SUBCRITICAL TEMPERATURE, TRANSVERSE PROPAGATOR */
    let mut plot = Plot2D::new();
    plot.set_xlabel("$p$ [GeV]");
    plot.set_ylabel("$\\Delta_{T}(p)$ [GeV$^{-2}$]");
    plot.set_path(
        &targetdir
            .join("gluon_propagator_t_low_t.png")
            .to_string_lossy(),
    );
    plot.set_domain(momenta.clone());

    plot.insert_image(vacuum_vals);
    let mut legends = vec![String::from("$T=0$ MeV")];

    for i in 0..=2 {
        let t = temperatures[i];
        let beta = 1. / t;
        let m = ms_t[i];
        let f0 = f0s_t[i];
        let norm = propagator_t_landau(om, renpoint, m, beta, f0).re * (renpoint * renpoint);
        let vals = momenta
            .iter()
            .map(|&p| {
                eprint!("Computing {p}... ");
                let res = propagator_t_landau(om, p, m, beta, f0).re / norm;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(vals);
        legends.push(format!("$T={}$ MeV", (t * 1000.) as u32));
    }
    plot.set_legend(legends.iter().map(|l| l.as_str()).collect());
    plot.savefig().expect("Could not save figure");

    /* SUPERCRITICAL TEMPERATURE, TRANSVERS PROPAGATOR */
    let mut plot = Plot2D::new();
    plot.set_xlabel("$p$ [GeV]");
    plot.set_ylabel("$\\Delta_{T}(p)$ [GeV$^{-2}$]");
    plot.set_path(
        &targetdir
            .join("gluon_propagator_t_high_t.png")
            .to_string_lossy(),
    );
    plot.set_domain(momenta.clone());
    let mut legends = vec![];
    for i in 3..=5 {
        let t = temperatures[i];
        let beta = 1. / t;
        let m = ms_t[i];
        let f0 = f0s_t[i];
        let norm = propagator_t_landau(om, renpoint, m, beta, f0).re * (renpoint * renpoint);
        let vals = momenta
            .iter()
            .map(|&p| {
                eprint!("Computing {p}... ");
                let res = propagator_t_landau(om, p, m, beta, f0).re / norm;
                eprintln!("Computed {res}.");
                res
            })
            .collect();
        plot.insert_image(vals);
        legends.push(format!("$T={}$ MeV", (t * 1000.) as u32));
    }
    plot.set_legend(legends.iter().map(|l| l.as_str()).collect());
    plot.savefig().expect("Could not save figure");
}
