use lazy_static::lazy_static;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use tc_long_prop::more_masses::*;

const D2W: usize = 2;

lazy_static! {
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("physical_tc");
}

#[derive(Clone)]
struct Config<'a> {
    // In units of m
    pbase: R,
    // In GeV
    om: R,
    // In GeV
    renpoint: R,
    // Adimensional
    f0: R,
    fieldconfig: &'a FieldConfig,
    correctedfieldconfig: Option<&'a FieldConfig>,
    label: &'a str,
    filename: &'a str,
    title: &'a str,
    temps: Vec<R>,
    chempots: Vec<R>,
    momenta: Vec<R>,
    renpoint2: R,
    correctedf0: Option<R>,
}

impl<'a> Config<'a> {
    fn new(
        temps: (R, R, R),
        chempots: (R, R, R),
        momenta: (R, R, R),
        om: R,
        renpoint: R,
        f0: R,
        fieldconfig: &'a FieldConfig,
        correctedfieldconfig: Option<&'a FieldConfig>,
        label: &'a str,
        filename: &'a str,
        title: &'a str,
    ) -> Self {
        let pbase = momenta.0;

        let temps = Self::compute_vec(temps);
        let chempots = Self::compute_vec(chempots);
        let momenta = Self::compute_vec(momenta);

        let renpoint2 = renpoint * renpoint;

        let correctedf0 = Self::compute_correctedf0(fieldconfig, correctedfieldconfig, f0);

        Self {
            pbase,
            om,
            renpoint,
            f0,
            fieldconfig,
            correctedfieldconfig,
            label,
            filename,
            title,
            temps,
            chempots,
            momenta,
            renpoint2,
            correctedf0,
        }
    }

    fn reset_temperatures(&mut self, temps: (R, R, R)) {
        self.temps = Self::compute_vec(temps);
    }

    fn reset_chemicalpotentials(&mut self, chpots: (R, R, R)) {
        self.chempots = Self::compute_vec(chpots);
    }

    fn reset_correctedfieldconfig(&mut self, correctedfieldconfig: Option<&'a FieldConfig>) {
        self.correctedfieldconfig = correctedfieldconfig;
        self.correctedf0 =
            Self::compute_correctedf0(self.fieldconfig, correctedfieldconfig, self.f0);
    }

    fn compute_vec(v: (R, R, R)) -> Vec<R> {
        let (vmin, vmax, dv) = v;
        let vrange = vmax - vmin;
        let n = (vrange / dv) as usize;
        (0..=n).map(|k| vmin + dv * (k as R)).collect()
    }

    fn compute_correctedf0<'b>(
        fieldconfig: &'b FieldConfig,
        correctedfieldconfig: Option<&'b FieldConfig>,
        f0: R,
    ) -> Option<R> {
        correctedfieldconfig.map(|cfc| {
            f0 + cfc
                .quarks
                .iter()
                .enumerate()
                .map(|(i, (nf, mqc))| (*nf as R) * (fieldconfig.quarks[i].1 / mqc).ln())
                .sum::<R>()
                * 4.
                / (9. * (fieldconfig.nc as R))
        })
    }
}

fn compute_propagators(config: &Config) {
    let fieldconfig = config.fieldconfig;

    let m = fieldconfig.gluon;
    let om = config.om;
    let (renpoint, renpoint2) = (config.renpoint, config.renpoint2);
    let f0 = config.f0;

    let (label, filename, title) = (config.label, config.filename, config.title);

    let renfac = renpoint2 / (m * m);
    let legends: Vec<String> = config
        .temps
        .iter()
        .map(|t| format!("$T/m={t:.4}$"))
        .collect();

    fs::create_dir_all(THIS_BASEDIR.as_path()).expect("Could not create base directory");

    for &mu in &config.chempots {
        let mut plot = Plot2D::new();
        plot.set_xlabel("$p/m$");
        plot.set_ylabel("$m^{2}\\,\\Delta_{L}(p)$");
        plot.set_title(&format!("{}, $\\mu={mu:.4}$ GeV", title));
        plot.set_domain(config.momenta.clone());
        plot.set_legend(legends.iter().map(|l| l.as_str()).collect());

        let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renfac;
        plot.insert_image(
            config
                .momenta
                .iter()
                .map(|p| {
                    let val = propagator_l_zero_temp_landau(om, p * m, mu, f0, fieldconfig).re / z;
                    eprintln!("Computed (T/m, mu, p/m) = (0.0000, {mu:.4}, {p:.4}) for {label}.");
                    val
                })
                .collect(),
        );

        config
            .temps
            .iter()
            .skip(1)
            .map(|t| {
                let beta = 1. / (t * m);
                let z =
                    propagator_l_landau(om, renpoint, beta, mu, f0, config.fieldconfig).re * renfac;
                config
                    .momenta
                    .par_iter()
                    .map(|p| {
                        let val =
                            propagator_l_landau(om, p * m, beta, mu, f0, config.fieldconfig).re / z;
                        eprintln!(
                            "Computed (T/m, mu, p/m) = ({t:.4}, {mu:.4}, {p:.4}) for {label}."
                        );
                        val
                    })
                    .collect()
            })
            .for_each(|vals: Vec<R>| {
                plot.insert_image(vals);
            });
        plot.set_path(&format!(
            "{}/{}_mu_{mu:.4}.png",
            THIS_BASEDIR.to_string_lossy(),
            filename
        ));
        plot.savefig().expect("Could not save figure");
    }
}

fn compute_ir_limit(config: &Config) {
    let fieldconfig = config.fieldconfig;

    let m = fieldconfig.gluon;
    let (pbase, om) = (config.pbase, config.om);
    let (renpoint, renpoint2) = (config.renpoint, config.renpoint2);
    let f0 = config.f0;

    let (label, filename, title) = (config.label, config.filename, config.title);

    fs::create_dir_all(THIS_BASEDIR.as_path()).expect("Could not create base directory");

    let mut plot = Plot2D::new();
    plot.set_domain(config.temps.clone());
    plot.set_xlabel("$T/m$");
    plot.set_ylabel("$m^{2}\\,\\Delta_{L}(p=0)$");
    plot.set_title(title);

    for &mu in &config.chempots {
        let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renpoint2;
        let mut vals =
            vec![propagator_l_zero_temp_landau(om, pbase * m, mu, f0, fieldconfig).re / z];
        eprintln!("Computed (T/m, mu) = (0.0000, {mu:.4}) for {label}.");
        vals.extend::<Vec<R>>(
            config
                .temps
                .par_iter()
                .skip(1)
                .map(|t| {
                    let beta = 1. / (t * m);
                    let z =
                        propagator_l_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renpoint2;
                    let val = propagator_l_landau(om, pbase * m, beta, mu, f0, fieldconfig).re / z;
                    eprintln!("Computed (T/m, mu) = ({t:.4}, {mu:.4}) for {label}.");
                    val
                })
                .collect(),
        );
        plot.insert_image(vals);
    }

    plot.set_path(&format!(
        "{}/ir_{filename}.png",
        THIS_BASEDIR.to_string_lossy()
    ));
    plot.savefig().expect("Could not save figure");
}

fn compute_phase_diagram(config: &Config) -> Vec<(R, R)> {
    let mut config = config.clone();

    let m = config.fieldconfig.gluon;
    let (pbase, om) = (config.pbase, config.om);
    let (renpoint, renpoint2) = (config.renpoint, config.renpoint2);

    let (label, filename, title) = (config.label, config.filename, config.title);

    let mut tcs = vec![];
    let mut correct = (true, config.correctedfieldconfig.is_some());

    for &mu in &config.chempots {
        let z = propagator_l_zero_temp_landau(om, renpoint, mu, config.f0, config.fieldconfig).re
            * renpoint2;
        let mut vals = vec![
            propagator_l_zero_temp_landau(om, pbase * m, mu, config.f0, config.fieldconfig).re / z,
        ];
        eprintln!("Computed (T/m, mu) = (0.0000, {mu:.4}) for {label}.");
        vals.extend::<Vec<R>>(
            config
                .temps
                .par_iter()
                .skip(1)
                .map(|t| {
                    let beta = 1. / (t * m);
                    let z =
                        propagator_l_landau(om, renpoint, beta, mu, config.f0, config.fieldconfig)
                            .re
                            * renpoint2;
                    let val =
                        propagator_l_landau(om, pbase * m, beta, mu, config.f0, config.fieldconfig)
                            .re
                            / z;
                    eprintln!("Computed (T/m, mu) = ({t:.4}, {mu:.4}) for {label}.");
                    val
                })
                .collect(),
        );
        let tc = find_critical_temperature(&config.temps, &vals);
        tcs.push(tc);
        // This is not super-precise in general. A better condition would
        // be to require tc != 0 but very small. Nonetheless, if there are
        // enough suddivisions in the mu variable, this is actually *more*
        // precise.
        if correct == (true, true) && tc == 0. {
            config.fieldconfig = config.correctedfieldconfig.unwrap();
            config.f0 = config.correctedf0.unwrap();
            correct.0 = false;
        }
    }

    let (cpl, d2w) = (config.chempots.len(), D2W);
    let dmu = config.chempots[1] - config.chempots[0];
    let mut tcsd2 = Vec::with_capacity(cpl - 2 * d2w);
    for i in d2w..cpl - d2w {
        let d2 = (tcs[i + d2w] + tcs[i - d2w] - 2. * tcs[i]) / (((d2w * d2w) as R) * dmu * dmu);
        tcsd2.push((config.chempots[i], d2));
    }

    fs::create_dir_all(THIS_BASEDIR.as_path()).expect("Could not create base directory");
    let mut outfile = BufWriter::new(
        fs::File::create(THIS_BASEDIR.join(format!("tc_{filename}.out")))
            .expect("Could not create output file"),
    );
    for (mu, tc) in config.chempots.iter().zip(&tcs) {
        writeln!(outfile, "{mu:.4}\t{tc:.4}").expect("Could not write value to file");
    }

    let mut outfile = BufWriter::new(
        fs::File::create(THIS_BASEDIR.join(format!("tcd2_{filename}.out")))
            .expect("Could not create output file"),
    );

    for (mu, tcd2) in &tcsd2 {
        writeln!(outfile, "{mu}\t{tcd2}").expect("Could not write value to file");
    }

    let mut plot = Plot2D::new();
    plot.set_domain(config.chempots.clone());
    plot.set_xlabel("$\\mu$ [GeV]");
    plot.set_ylabel("$T_{c}/m$");
    plot.set_title(title);
    plot.insert_image(tcs.clone());
    plot.set_path(&format!(
        "{}/tc_{filename}.png",
        THIS_BASEDIR.to_string_lossy()
    ));
    plot.savefig().expect("Could not save figure");

    let mut plot = Plot2D::new();
    plot.set_domain(tcsd2.iter().map(|(mu, _)| *mu).collect());
    plot.set_xlabel("$\\mu$ [GeV]");
    plot.set_ylabel("$\\partial^{2}(T_{c}/m)/\\partial \\mu^{2}$ [GeV$^{-2}$]");
    plot.set_ylim((-10., 10.));
    plot.set_title(title);
    plot.insert_image(tcsd2.iter().map(|(_, t)| *t).collect());
    plot.set_path(&format!(
        "{}/tcd2_{filename}.png",
        THIS_BASEDIR.to_string_lossy()
    ));
    plot.savefig().expect("Could not save figure");

    config
        .chempots
        .iter()
        .zip(tcs.iter())
        .map(|(&mu, &tc)| (mu, tc))
        .collect()
}

fn main() {
    /* CRITICAL TEMPERATURE AS A FUNCTION OF THE CHEMICAL POTENTIAL
    FOR A PLAUSIBLE PHYSICAL CONFIGURATION, FIRST WITH NF = 2 + 1,
    THEN WITH NF = 2 + 1 + 1 */

    // In GeV
    let m = 0.656;
    let om = 1E-5;
    let renpoint = 4.;
    let (mq1, mq2, mq3) = (0.350, 0.450, 1.5);
    let (mq1c, mq2c, mq3c) = (0.003, 0.090, 1.2);

    // Adimensional
    let f0 = -0.876;

    // In units of m
    let fewtemps = (0., 0.25, 0.05);
    let moretemps = (0., 0.25, 0.01);
    let manytemps = (0., 0.25, 0.0001);

    // In GeV
    let fewchempots = (0., 1., 0.1);
    let morechempots = (0., 1., 0.025);
    let manychempots = (0., 1., 0.01);

    // In units of m
    let momenta = (0.01, 3., 0.01);

    /* NF = 2 + 1 */
    let fieldconfig = FieldConfig::new(3, m, vec![(2, mq1), (1, mq2)]);
    let correctedfieldconfig = FieldConfig::new(3, m, vec![(2, mq1c), (1, mq2c)]);
    let mut config = Config::new(
        fewtemps,
        fewchempots,
        momenta,
        om,
        renpoint,
        f0,
        &fieldconfig,
        None,
        "nf = 2+1",
        "nf_2+1",
        "$n_{{F}}=2+1$",
    );
    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 ***");
    compute_propagators(&config);
    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 AROUND THE QUARK MASSES ***");
    let oldfilename = config.filename;
    config.reset_temperatures((0., 0.05, 0.01));
    config.reset_chemicalpotentials((mq1 - 0.05, mq1 + 0.05, 0.01));
    config.filename = "nf_2+1_m1";
    compute_propagators(&config);
    config.reset_temperatures((0., 0.1, 0.02));
    config.reset_chemicalpotentials((mq2 - 0.05, mq2 + 0.05, 0.01));
    config.filename = "nf_2+1_m2";
    compute_propagators(&config);
    config.filename = oldfilename;
    eprintln!("*** COMPUTING PROPAGATORS's IR LIMIT AT NF = 2 + 1 ***");
    config.reset_temperatures(moretemps);
    config.reset_chemicalpotentials(morechempots);
    compute_ir_limit(&config);
    config.reset_temperatures(manytemps);
    config.reset_chemicalpotentials(manychempots);
    eprintln!("*** COMPUTING PHASE DIAGRAM AT NF = 2 + 1 ***");
    compute_phase_diagram(&config);
    eprintln!("*** COMPUTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 ***");
    config.reset_correctedfieldconfig(Some(&correctedfieldconfig));
    config.filename = "nf_2+1_corrected";
    let cpd = compute_phase_diagram(&config);
    let pdf = parametrize_phase_diagram(
        &cpd,
        8,
        THIS_BASEDIR
            .join(format!("tc_{}_fit.png", config.filename))
            .to_str(),
    );
    println!("{pdf:?}"); //DEBUG

    /* NF = 2 + 1 + 1 */
    let fieldconfig = FieldConfig::new(3, m, vec![(2, mq1), (1, mq2), (1, mq3)]);
    let correctedfieldconfig = FieldConfig::new(3, m, vec![(2, mq1c), (1, mq2c), (1, mq3c)]);
    let mut config = Config::new(
        fewtemps,
        fewchempots,
        momenta,
        om,
        renpoint,
        f0,
        &fieldconfig,
        None,
        "nf = 2+1+1",
        "nf_2+1+1",
        "$n_{{F}}=2+1+1$",
    );
    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 + 1 ***");
    compute_propagators(&config);
    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 AROUND THE QUARK MASSES ***");
    let oldfilename = config.filename;
    config.reset_temperatures((0., 0.05, 0.01));
    config.reset_chemicalpotentials((mq1 - 0.05, mq1 + 0.05, 0.01));
    config.filename = "nf_2+1+1_m1";
    compute_propagators(&config);
    config.reset_temperatures((0., 0.1, 0.02));
    config.reset_chemicalpotentials((mq2 - 0.05, mq2 + 0.05, 0.01));
    config.filename = "nf_2+1+1_m2";
    compute_propagators(&config);
    config.reset_temperatures((0., 0.2, 0.04));
    config.reset_chemicalpotentials((mq3 - 0.3, mq3 + 1.2, 0.15));
    config.filename = "nf_2+1+1_m3";
    compute_propagators(&config);
    config.filename = oldfilename;
    eprintln!("*** COMPUTING PROPAGATORS's IR LIMIT AT NF = 2 + 1 ***");
    config.reset_temperatures(moretemps);
    config.reset_chemicalpotentials(morechempots);
    compute_ir_limit(&config);
    config.reset_temperatures(manytemps);
    config.reset_chemicalpotentials(manychempots);
    eprintln!("*** COMPUTING PHASE DIAGRAM AT NF = 2 + 1 + 1 ***");
    config.reset_chemicalpotentials((0., 3., 0.01));
    compute_phase_diagram(&config);
    eprintln!("*** COMPUTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 + 1 ***");
    config.reset_correctedfieldconfig(Some(&correctedfieldconfig));
    config.filename = "nf_2+1+1_corrected";
    compute_phase_diagram(&config);
}
