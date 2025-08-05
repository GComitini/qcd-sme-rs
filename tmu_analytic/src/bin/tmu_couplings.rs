use qcd_sme::qcd::thermal::gluon::{
    propagator_t_landau_w_field_config as qcd_propagator,
    propagator_t_zero_temp_landau_w_field_config as qcd_propagator_zero_temp,
};
use qcd_sme::{qcd::FieldConfig, types::NCTYPE, R};
use std::f64::consts::PI;

const MG: R = 0.656;
const NC: NCTYPE = 3;
const F0: R = -0.876;
const PREN: R = 4.;
const OMEPS: R = 1E-3;

fn alpha_s_from_f0_and_pren(f0: R, t: R, pren: R, config: &FieldConfig) -> R {
    4. * PI / (3. * (NC as R))
        * pren
        * pren
        * if t == 0. {
            qcd_propagator_zero_temp(OMEPS, pren, 0., f0, config).re
        } else {
            qcd_propagator(OMEPS, pren, 1. / t, 0., f0, config).re
        }
}

fn alpha_s_from_f0(f0: R, t: R, config: &FieldConfig) -> R {
    alpha_s_from_f0_and_pren(f0, t, PREN, config)
}

fn f0_from_alpha_s_and_pren(a_s: R, t: R, pren: R, config: &FieldConfig) -> R {
    4. * PI / (3. * (NC as R) * a_s)
        - 1. / (pren
            * pren
            * if t == 0. {
                qcd_propagator_zero_temp(OMEPS, pren, 0., 0., config).re
            } else {
                qcd_propagator(OMEPS, pren, 1. / t, 0., 0., config).re
            })
}

fn f0_from_alpha_s(a_s: R, t: R, config: &FieldConfig) -> R {
    f0_from_alpha_s_and_pren(a_s, t, PREN, config)
}

fn main() {
    /* I. Pure Yang-Mills theory */
    // IA. Massive
    let config = FieldConfig::new(NC, MG, vec![]);
    println!("* Pure Yang-Mills theory");
    println!(
        "T = 0, F0 = {F0} => a_s = {}",
        alpha_s_from_f0(F0, 0., &config)
    );
    let a_s = 0.5167177565590644;
    println!(
        "T = 0, a_s = {a_s} => F0 = {}",
        f0_from_alpha_s(a_s, 0., &config)
    );
    println!("");

    // IB. Massless
    qcd_sme::consts::set_default_max_iter_integral(20);
    let config = FieldConfig::new(NC, 0.005, vec![]);
    let (t, f0ml, pren) = (0.1, 25., 0.5);
    println!("* Massless pure Yang-Mills theory at p_ren = {pren} GeV");
    println!(
        "T = {t}, F0 = {f0ml} => a_s = {}",
        alpha_s_from_f0_and_pren(f0ml, t, pren, &config)
    );
    let a_s = 0.04153253103565679;
    println!(
        "T = {t}, a_s = {a_s} => F0 = {}",
        f0_from_alpha_s_and_pren(a_s, t, pren, &config)
    );
    println!("");
    qcd_sme::consts::set_default_max_iter_integral(50);

    /* II. Full QCD */
    // IIA. Pure Yang-Mills parameters
    let config = FieldConfig::new(NC, MG, vec![(2, 0.35), (1, 0.45)]);
    let f00 = F0
        + config
            .quarks
            .iter()
            .map(|(nf, mq)| (*nf as R) * (config.gluon / mq).ln())
            .sum::<R>()
            * 4.
            / (9. * (config.nc as R));
    println!("* Full QCD, quarkconfig = [(2, 0.35), (1, 0.45)], native F0 = {F0}");
    println!(
        "T = 0, F0 = {f00} => a_s = {}",
        alpha_s_from_f0(f00, 0., &config)
    );
    let a_s = 0.618573343225959;
    println!(
        "T = 0, a_s = {a_s} => F0 = {}",
        f0_from_alpha_s(a_s, 0., &config)
    );
    println!("Note: on the last two lines, F0 is given in qcd_sme library conventions");
    println!("");

    // IIB. Nf = 2 lowest temperature fit parameters @ mq = 0.4 GeV
    const F0QCD: R = -0.506;
    let config = FieldConfig::new(NC, MG, vec![(2, 0.4)]);
    // F0QCD is given in library convention. Compute the native f0
    let f0n = F0QCD
        - config
            .quarks
            .iter()
            .map(|(nf, mq)| (*nf as R) * (config.gluon / mq).ln())
            .sum::<R>()
            * 4.
            / (9. * (config.nc as R));
    println!("* Full QCD, quarkconfig = [(2, 0.4)], native F0 = {f0n}");
    println!(
        "T = 0, F0 = {F0QCD} => a_s = {}",
        alpha_s_from_f0(F0QCD, 0., &config)
    );
    let a_s = 0.5312208877232416;
    println!(
        "T = 0, a_s = {a_s} => F0 = {}",
        f0_from_alpha_s(a_s, 0., &config)
    );
    println!("Note: on the last two lines, F0 is given in qcd_sme library conventions");
}
