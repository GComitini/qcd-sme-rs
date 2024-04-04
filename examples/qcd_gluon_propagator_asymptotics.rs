use peroxide::fuga::*;
use qcd_sme::{
    consts::{
        get_number_of_colors, get_number_of_fermions, set_default_tol_integral,
        set_number_of_colors, set_number_of_fermions,
    },
    qcd::thermal::gluon::propagator_l_landau,
    R,
};
use std::f64::consts::PI;

// Logarithm of the HTL mass squared divided by the coupling constant squared
// modulo a constant prefactor
fn asymptotics_known(t: R, mu: R, nc: u32, nf: u32, rem: R) -> R {
    let nf = nf as R;
    let nc = nc as R;
    ((nc + nf / 2.) * t * t + 3. * nf * mu * mu / (2. * PI * PI)).ln() + rem
}

// For the love of god, do not use this function while multithreading
fn asymptotics_test(t: R, mu: R, nc: u32, nf: u32) -> R {
    let (nc_old, nf_old) = (get_number_of_colors(), get_number_of_fermions());
    set_number_of_colors(nc);
    set_number_of_fermions(nf);
    let (om, p, m, beta, f0) = (0.001, 0.001, 0.3, 1. / t, 1.5);
    // Our normalization is nc-dependent: remove this dependence by dividing by nc.
    let res = -(propagator_l_landau(om, p, m, beta, mu, f0).re / nc as R).ln();
    set_number_of_colors(nc_old);
    set_number_of_fermions(nf_old);
    res
}

/* In the limit of large temperatures and large chemical potentials (with respect to
  the momenta and masses), the gluon propagator approaches Z/(p^2 + mg^2), where mg^2
  is the operand of the logarithm in "asymptotics_known" (modulo a factor of g^2/3).
  If our definitions are correct, for large temperatures and, separately, for large
  chemical potentials, "asymptotics_known" and "asymptotics_test" will differ by a
  constant which is asymptotically independent from nc, nf, t and mu. This constant
  can be tuned using the "rem" argument of asymptotics_known. If tuning "rem" is
  enough to make "asymptotics_known" and "asymptotics_test" match, then our definitions
  are asymptotically correct.
*/
fn main() {
    set_default_tol_integral(1E-4);

    println!("*** TEMPERATURE COMPARISONS ***");

    let mut plot = peroxide::util::plot::Plot2D::new();
    let temps = vec![8., 12., 15., 20.];
    let mu = 3.;
    plot.set_domain(temps.clone());
    plot.set_path("target/qcd_gluon_propagator_asymptotics_temperature.png");
    plot.set_xlabel("$$T$$");
    let mut legends = Vec::new();

    for nf in [1, 2] {
        for nc in [2, 3, 5] {
            let (mut y1, mut y2) = (vec![], vec![]);
            for t in &temps {
                print!("Computing nf = {nf}, nc = {nc}, T = {t}... ");
                y1.push(asymptotics_known(*t, mu, nc, nf, 2.02));
                y2.push(asymptotics_test(*t, mu, nc, nf));
                println!("First diff: {}.", y2[0] - y1[0]);
            }
            plot.insert_image(y1);
            legends.push(format!("nf = {nf}, nc = {nc} (known)"));
            plot.insert_image(y2);
            legends.push(format!("nf = {nf}, nc = {nc} (test)"));
        }
    }
    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.savefig().expect("WTF");

    println!("\n*** CHEMICAL POTENTIAL COMPARISONS ***");

    let mut plot = peroxide::util::plot::Plot2D::new();
    let chempots = vec![8., 12., 15., 20.];
    let t = 3.;
    plot.set_domain(chempots.clone());
    plot.set_path("target/qcd_gluon_propagator_asymptotics_chempot.png");
    plot.set_xlabel("$$\\mu$$");
    let mut legends = Vec::new();

    for nf in [1, 2] {
        for nc in [2, 3, 5] {
            let (mut y1, mut y2) = (vec![], vec![]);
            for mu in &chempots {
                print!("Computing nf = {nf}, nc = {nc}, mu = {mu}... ");
                y1.push(asymptotics_known(t, *mu, nc, nf, 2.02));
                y2.push(asymptotics_test(t, *mu, nc, nf));
                println!("First diff: {}.", y2[0] - y1[0]);
            }
            plot.insert_image(y1);
            legends.push(format!("nf = {nf}, nc = {nc} (known)"));
            plot.insert_image(y2);
            legends.push(format!("nf = {nf}, nc = {nc} (test)"));
        }
    }
    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.savefig().expect("WTF");
}
