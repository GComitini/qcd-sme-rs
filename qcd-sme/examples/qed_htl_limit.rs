use peroxide::fuga::Plot;
use peroxide::util::plot::Plot2D;
use qcd_sme::low_level::oneloop::thermal::gluon::{
    polarization_quark_l_thermal_part_landau, polarization_quark_t_thermal_part_landau,
};
use qcd_sme::{Num, C, I, R};
use rayon::prelude::*;
use std::fs;

const ME: R = 0.001;

// Longitudinal QED polarization in the HTL limit divided by e^2 (positive)
fn f_htl<T: Num>(om: T, p: R, t: R) -> C {
    let k2 = -om * om - p * p;
    let p2 = p * p;
    let k0 = om * I;
    -k2 * (t * t / 3.) / p2 * (1. - k0 / (2. * p) * ((k0 + p) / (k0 - p)).ln())
}

// Transverse QED polarization in the HTL limit divided by e^2 (positive)
fn g_htl<T: Num>(om: T, p: R, t: R) -> C {
    t * t / 6. - f_htl(om, p, t) / 2.
}

// Longitudinal QED polarization as computed from an analogous term in QCD,
// divided by e^2. 6 = 2 x 3, where 3 = Nc and 2 comes from g^2 -> 2 e^2
fn f_ours<T: Num>(om: T, p: R, t: R) -> C {
    -6. * polarization_quark_l_thermal_part_landau(om, p, ME, 1. / t, 0.)
}

// Transverse QED polarization as computed from an analogous term in QCD,
// divided by e^2. 6 = 2 x 3, where 3 = Nc and 2 comes from g^2 -> 2 e^2
fn g_ours<T: Num>(om: T, p: R, t: R) -> C {
    -6. * polarization_quark_t_thermal_part_landau(om, p, ME, 1. / t, 0.)
}

fn main() {
    let om = 0.01;
    let ts = [1., 2., 5., 10.];

    let (pmin, pmax) = (0.01, 0.5);
    let dp = 0.02;
    let prange = pmax - pmin;
    let n = (prange / dp) as usize;
    let momenta: Vec<R> = (0..=n).map(|k| pmin + dp * (k as R)).collect();

    fs::create_dir_all("target/qed_htl_limit").expect("could not create plots directory");

    for t in ts {
        let mut legends: Vec<String> = vec![];

        let mut plot_re = Plot2D::new();
        plot_re.set_xlabel("$p$ [MeV]");
        plot_re.set_ylabel("Re{$\\Delta(p)$} [MeV$^{-2}$]");
        plot_re.set_title(&format!("$\\omega={om}$ MeV, $m_{{e}}={ME}$ MeV"));
        plot_re.set_domain(momenta.clone());
        plot_re.set_path(&format!("target/qed_htl_limit/T{:.2}_re.png", t));

        let mut plot_im = Plot2D::new();
        plot_im.set_xlabel("$p$ [MeV]");
        plot_im.set_ylabel("Im{$\\Delta(p)$} [MeV$^{-2}$]");
        plot_im.set_title(&format!("$\\omega={om}$ MeV, $m_{{e}}={ME}$ MeV"));
        plot_im.set_domain(momenta.clone());
        plot_im.set_path(&format!("target/qed_htl_limit/T{:.2}_im.png", t));

        let vals: Vec<C> = momenta
            .iter()
            .map(|&p| {
                eprintln!("Computing long. HTL (T,p) = ({t},{p})...");
                f_htl(om, p, t)
            })
            .collect();
        plot_re.insert_image(vals.iter().map(|v| v.re).collect());
        plot_im.insert_image(vals.iter().map(|v| v.im).collect());
        legends.push(format!("$T={t}$ MeV (long., HTL)"));

        let vals: Vec<C> = momenta
            .par_iter()
            .map(|&p| {
                eprintln!("Computing long. ours (T,p) = ({t},{p})...");
                f_ours(om, p, t)
            })
            .collect();
        plot_re.insert_image(vals.iter().map(|v| v.re).collect());
        plot_im.insert_image(vals.iter().map(|v| v.im).collect());
        legends.push(format!("$T={t}$ MeV (long., ours)"));

        let vals: Vec<C> = momenta
            .iter()
            .map(|&p| {
                eprintln!("Computing trans. HTL (T,p) = ({t},{p})...");
                g_htl(om, p, t)
            })
            .collect();
        plot_re.insert_image(vals.iter().map(|v| v.re).collect());
        plot_im.insert_image(vals.iter().map(|v| v.im).collect());
        legends.push(format!("$T={t}$ MeV (trans., HTL)"));

        let vals: Vec<C> = momenta
            .par_iter()
            .map(|&p| {
                eprintln!("Computing trans. ours (T,p) = ({t},{p})...");
                g_ours(om, p, t)
            })
            .collect();
        plot_re.insert_image(vals.iter().map(|v| v.re).collect());
        plot_im.insert_image(vals.iter().map(|v| v.im).collect());
        legends.push(format!("$T={t}$ MeV (trans., ours)"));

        plot_re.set_legend(legends.iter().map(|l| l.as_str()).collect());
        plot_re.savefig().expect("Could not save figure");

        plot_im.set_legend(legends.iter().map(|l| l.as_str()).collect());
        plot_im.savefig().expect("Could not save figure");
    }
}
