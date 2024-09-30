use config::Config;
use lazy_static::lazy_static;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use tc_long_prop::more_masses::*;

// Width (in points) for numerical second derivative
const D2W: usize = 2;

lazy_static! {
    // Base directory for saving output of this binary
    static ref THIS_BASEDIR: PathBuf = PathBuf::from(BASEDIR).join("physical_tc");
}

mod config {
    use crate::*;

    // Public fields can be modified independent of other fields so
    // they do not require special setter methods. We don't currently
    // reset self.momenta so we don't need to make pbase private
    #[derive(Clone)]
    pub struct Config<'a> {
        // In units of m
        pub pbase: R,
        // In GeV
        pub om: R,
        pub fieldconfig: &'a FieldConfig,
        pub correctedfieldconfig: Option<&'a FieldConfig>,
        pub label: &'a str,
        pub filename: &'a str,
        pub title: &'a str,
        // In GeV
        renpoint: R,
        // Adimensional
        f0: R,
        temps: Vec<R>,
        chempots: Vec<R>,
        momenta: Vec<R>,
        renpoint2: R,
        correctedf0: Option<R>,
        phase_boundary: Option<&'a [(R, R)]>,
    }

    impl<'a> Config<'a> {
        pub fn new(
            // (min, max, delta)
            temps: (R, R, R),
            chempots: (R, R, R),
            momenta: (R, R, R),
            om: R,
            renpoint: R,
            f0: R,
            fieldconfig: &'a FieldConfig,
            label: &'a str,
            filename: &'a str,
            title: &'a str,
        ) -> Self {
            let pbase = momenta.0;

            let temps = Self::compute_vec(temps);
            let chempots = Self::compute_vec(chempots);
            let momenta = Self::compute_vec(momenta);

            let renpoint2 = renpoint * renpoint;

            Self {
                pbase,
                om,
                fieldconfig,
                correctedfieldconfig: None,
                label,
                filename,
                title,
                renpoint,
                f0,
                temps,
                chempots,
                momenta,
                renpoint2,
                correctedf0: None,
                phase_boundary: None,
            }
        }

        pub fn renpoint(&self) -> R {
            self.renpoint
        }

        pub fn renpoint2(&self) -> R {
            self.renpoint2
        }

        pub fn temps(&self) -> &Vec<R> {
            &self.temps
        }

        pub fn chempots(&self) -> &Vec<R> {
            &self.chempots
        }

        pub fn momenta(&self) -> &Vec<R> {
            &self.momenta
        }

        pub fn phase_boundary(&self) -> Option<&[(R, R)]> {
            self.phase_boundary
        }

        pub fn reset_temperatures(&mut self, temps: (R, R, R)) {
            self.temps = Self::compute_vec(temps);
        }

        pub fn reset_chemicalpotentials(&mut self, chempots: (R, R, R)) {
            self.chempots = Self::compute_vec(chempots);
        }

        pub fn reset_correctedfieldconfig(
            &mut self,
            correctedfieldconfig: Option<&'a FieldConfig>,
        ) {
            self.correctedfieldconfig = correctedfieldconfig;
            self.set_correctedf0();
        }

        pub fn reset_phase_boundary(&mut self, phase_boundary: Option<&'a [(R, R)]>) {
            // If any, reduce the phase boundary before storing it
            self.phase_boundary = phase_boundary.map(reduce_phase_boundary);
        }

        pub fn maybe_corrected_data(&self, mu: R, t: R) -> (&'a FieldConfig, R) {
            match (
                self.correctedfieldconfig,
                self.correctedf0,
                self.phase_boundary,
            ) {
                (Some(cfc), Some(cf0), Some(pb)) => {
                    if is_deconfined_phase(mu, t, pb) {
                        (cfc, cf0)
                    } else {
                        (self.fieldconfig, self.f0)
                    }
                }
                (None, None, None) => (self.fieldconfig, self.f0),
                _ => panic!("Inconsistent 'corrected' configuration"),
            }
        }

        fn set_correctedf0(&mut self) {
            // Corrects f0: adds it appropriate ln(mq/mqc) terms so that the vacuum
            // quark loop does not diverge in the mq -> 0 limit due to a bad choice
            // of renormalization constants.
            self.correctedf0 = self.correctedfieldconfig.map(|cfc| {
                self.f0
                    + cfc
                        .quarks
                        .iter()
                        .enumerate()
                        .map(|(i, (nf, mqc))| {
                            let mq = self.fieldconfig.quarks[i].1;
                            (*nf as R) * (mq / mqc).ln()
                        })
                        .sum::<R>()
                        * 4.
                        / (9. * (self.fieldconfig.nc as R))
            });
        }

        fn compute_vec(v: (R, R, R)) -> Vec<R> {
            let (vmin, vmax, dv) = v;
            let vrange = vmax - vmin;
            let n = (vrange / dv).round() as usize;
            (0..=n).map(|k| vmin + dv * (k as R)).collect()
        }
    }
}

fn compute_propagators(config: &Config) {
    let m = config.fieldconfig.gluon;
    let om = config.om;
    let (renpoint, renpoint2) = (config.renpoint(), config.renpoint2());
    let (label, filename, title) = (config.label, config.filename, config.title);

    let renfac = renpoint2 / (m * m);
    let legends: Vec<String> = config
        .temps()
        .iter()
        .map(|t| format!("$T/m={t:.4}$"))
        .collect();

    fs::create_dir_all(THIS_BASEDIR.join("data")).expect("Could not create base directory");

    for &mu in config.chempots() {
        let mut plot = Plot2D::new();
        plot.set_xlabel("$p/m$");
        plot.set_ylabel("$m^{2}\\,\\Delta_{L}(p)$");
        plot.set_title(&format!("{}, $\\mu={mu:.4}$ GeV", title));
        plot.set_domain(config.momenta().clone());
        plot.set_legend(legends.iter().map(|l| l.as_str()).collect());

        // If we want the corrected propagator (i.e. a propagator with different
        // quark masses, to mimick chiral symmetry restoration) we reset the field
        // configuration and f0 depending on the chemical potential and temperature.
        // We must (potentially) reset it for each (mu, T) pair before computing the
        // renormalization factors z
        let (fieldconfig, f0) = config.maybe_corrected_data(mu, 0.);
        let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renfac;
        let propvals: Vec<f64> = config
            .momenta()
            .iter()
            .map(|p| {
                let val = propagator_l_zero_temp_landau(om, p * m, mu, f0, fieldconfig).re / z;
                eprintln!("Computed (T/m, mu, p/m) = (0.0000, {mu:.4}, {p:.4}) for {label}.");
                val
            })
            .collect();
        let outfilename = format!(
            "{}_glulongprop_t_0.000_mu_{mu:.4}.out",
            THIS_BASEDIR.join("data").join(filename).to_string_lossy()
        );
        let mut outfile = BufWriter::new(
            fs::File::create(&outfilename)
                .unwrap_or_else(|_| panic!("could not create {outfilename}")),
        );
        config
            .momenta()
            .iter()
            .zip(propvals.iter())
            .for_each(|(p, d)| {
                writeln!(outfile, "{p}\t{d}")
                    .unwrap_or_else(|_| panic!("could not write to {outfilename}"));
            });
        plot.insert_image(propvals);

        config
            .temps()
            .iter()
            .skip(1)
            .map(|t| {
                let beta = 1. / (t * m);
                let (fieldconfig, f0) = config.maybe_corrected_data(mu, *t);
                let z = propagator_l_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renfac;
                let outfilename = format!(
                    "{}_glulongprop_t_{:.3}_mu_{mu:.4}.out",
                    THIS_BASEDIR.join("data").join(filename).to_string_lossy(),
                    t
                );
                let mut outfile = BufWriter::new(
                    fs::File::create(&outfilename)
                        .unwrap_or_else(|_| panic!("could not create {outfilename}")),
                );
                let propvals: Vec<f64> = config
                    .momenta()
                    .par_iter()
                    .map(|p| {
                        let val = propagator_l_landau(om, p * m, beta, mu, f0, fieldconfig).re / z;
                        eprintln!(
                            "Computed (T/m, mu, p/m) = ({t:.4}, {mu:.4}, {p:.4}) for {label}."
                        );
                        val
                    })
                    .collect();
                config
                    .momenta()
                    .iter()
                    .zip(propvals.iter())
                    .for_each(|(p, d)| {
                        writeln!(outfile, "{p}\t{d}")
                            .unwrap_or_else(|_| panic!("could not write to {outfilename}"));
                    });
                propvals
            })
            .for_each(|vals: Vec<R>| {
                plot.insert_image(vals);
            });
        plot.set_path(&format!(
            "{}_mu_{mu:.4}.png",
            THIS_BASEDIR.join(filename).to_string_lossy()
        ));
        plot.savefig().expect("Could not save figure");
    }
}

fn compute_transverse_propagators(config: &Config) {
    let m = config.fieldconfig.gluon;
    let om = config.om;
    let (renpoint, renpoint2) = (config.renpoint(), config.renpoint2());
    let (label, filename, title) = (config.label, config.filename, config.title);

    let renfac = renpoint2 / (m * m);
    let legends: Vec<String> = config
        .temps()
        .iter()
        .map(|t| format!("$T/m={t:.4}$"))
        .collect();

    fs::create_dir_all(THIS_BASEDIR.join("data")).expect("Could not create base directory");

    for &mu in config.chempots() {
        let mut plot = Plot2D::new();
        plot.set_xlabel("$p/m$");
        plot.set_ylabel("$m^{2}\\,\\Delta_{T}(p)$");
        plot.set_title(&format!("{}, $\\mu={mu:.4}$ GeV", title));
        plot.set_domain(config.momenta().clone());
        plot.set_legend(legends.iter().map(|l| l.as_str()).collect());

        // If we want the corrected propagator (i.e. a propagator with different
        // quark masses, to mimick chiral symmetry restoration) we reset the field
        // configuration and f0 depending on the chemical potential and temperature.
        // We must (potentially) reset it for each (mu, T) pair before computing the
        // renormalization factors z
        let (fieldconfig, f0) = config.maybe_corrected_data(mu, 0.);
        let z = propagator_t_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renfac;
        let propvals: Vec<f64> = config
            .momenta()
            .iter()
            .map(|p| {
                let val = propagator_t_zero_temp_landau(om, p * m, mu, f0, fieldconfig).re / z;
                eprintln!("Computed (T/m, mu, p/m) = (0.0000, {mu:.4}, {p:.4}) for {label}.");
                val
            })
            .collect();
        let outfilename = format!(
            "{}_glutransprop_t_0.000_mu_{mu:.4}.out",
            THIS_BASEDIR.join("data").join(filename).to_string_lossy()
        );
        let mut outfile = BufWriter::new(
            fs::File::create(&outfilename)
                .unwrap_or_else(|_| panic!("could not create {outfilename}")),
        );
        config
            .momenta()
            .iter()
            .zip(propvals.iter())
            .for_each(|(p, d)| {
                writeln!(outfile, "{p}\t{d}")
                    .unwrap_or_else(|_| panic!("could not write to {outfilename}"));
            });
        plot.insert_image(propvals);

        config
            .temps()
            .iter()
            .skip(1)
            .map(|t| {
                let beta = 1. / (t * m);
                let (fieldconfig, f0) = config.maybe_corrected_data(mu, *t);
                let z = propagator_t_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renfac;
                let outfilename = format!(
                    "{}_glutransprop_t_{:.3}_mu_{mu:.4}.out",
                    THIS_BASEDIR.join("data").join(filename).to_string_lossy(),
                    t
                );
                let mut outfile = BufWriter::new(
                    fs::File::create(&outfilename)
                        .unwrap_or_else(|_| panic!("could not create {outfilename}")),
                );
                let propvals: Vec<f64> =config
                    .momenta()
                    .par_iter()
                    .map(|p| {
                        let val = propagator_t_landau(om, p * m, beta, mu, f0, fieldconfig).re / z;
                        eprintln!(
                            "Computed (T/m, mu, p/m) = ({t:.4}, {mu:.4}, {p:.4}) for {label} (transverse)."
                        );
                        val
                    })
                    .collect();
                config
                    .momenta()
                    .iter()
                    .zip(propvals.iter())
                    .for_each(|(p, d)| {
                        writeln!(outfile, "{p}\t{d}")
                            .unwrap_or_else(|_| panic!("could not write to {outfilename}"));
                    });
                propvals
            })
            .for_each(|vals: Vec<R>| {
                plot.insert_image(vals);
            });
        plot.set_path(&format!(
            "{}_transverse_mu_{mu:.4}.png",
            THIS_BASEDIR.join(filename).to_string_lossy()
        ));
        plot.savefig().expect("Could not save figure");
    }
}

fn compute_ir_limit(config: &Config) {
    let m = config.fieldconfig.gluon;
    let (pbase, om) = (config.pbase, config.om);
    let (renpoint, renpoint2) = (config.renpoint(), config.renpoint2());
    let (label, filename, title) = (config.label, config.filename, config.title);

    let renfac = renpoint2 / (m * m);

    fs::create_dir_all(THIS_BASEDIR.as_path()).expect("Could not create base directory");

    let mut plot = Plot2D::new();
    plot.set_domain(config.temps().clone());
    plot.set_xlabel("$T/m$");
    plot.set_ylabel("$m^{2}\\,\\Delta_{L}(p=0)$");
    plot.set_title(title);

    for &mu in config.chempots() {
        let (fieldconfig, f0) = config.maybe_corrected_data(mu, 0.);
        let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renfac;
        let mut vals =
            vec![propagator_l_zero_temp_landau(om, pbase * m, mu, f0, fieldconfig).re / z];
        eprintln!("Computed (T/m, mu) = (0.0000, {mu:.4}) for {label}.");

        vals.extend::<Vec<R>>(
            config
                .temps()
                .par_iter()
                .skip(1)
                .map(|t| {
                    let beta = 1. / (t * m);
                    let (fieldconfig, f0) = config.maybe_corrected_data(mu, 0.);
                    let z =
                        propagator_l_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renfac;
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
    let m = config.fieldconfig.gluon;
    let (pbase, om) = (config.pbase, config.om);
    let (renpoint, renpoint2) = (config.renpoint(), config.renpoint2());
    let (label, filename, title) = (config.label, config.filename, config.title);

    let renfac = renpoint2 / (m * m);

    let mut tcs = vec![];

    for &mu in config.chempots() {
        let (fieldconfig, f0) = config.maybe_corrected_data(mu, 0.);
        let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renfac;
        let mut vals =
            vec![propagator_l_zero_temp_landau(om, pbase * m, mu, f0, fieldconfig).re / z];
        eprintln!("Computed (T/m, mu) = (0.0000, {mu:.4}) for {label}.");

        vals.extend::<Vec<R>>(
            config
                .temps()
                .par_iter()
                .skip(1)
                .map(|t| {
                    let beta = 1. / (t * m);
                    let (fieldconfig, f0) = config.maybe_corrected_data(mu, *t);
                    let z =
                        propagator_l_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renfac;
                    let val = propagator_l_landau(om, pbase * m, beta, mu, f0, fieldconfig).re / z;
                    eprintln!("Computed (T/m, mu) = ({t:.4}, {mu:.4}) for {label}.");
                    val
                })
                .collect(),
        );
        let tc = find_critical_temperature(config.temps(), &vals);
        tcs.push(tc);
    }

    let (cpl, d2w) = (config.chempots().len(), D2W);
    let dmu = config.chempots()[1] - config.chempots()[0];
    let mut tcsd2 = Vec::with_capacity(cpl - 2 * d2w);
    let dmu2r = ((d2w * d2w) as R) * dmu * dmu;
    for i in d2w..cpl - d2w {
        let d2 = (tcs[i + d2w] + tcs[i - d2w] - 2. * tcs[i]) / dmu2r;
        tcsd2.push((config.chempots()[i], d2));
    }

    fs::create_dir_all(THIS_BASEDIR.as_path()).expect("Could not create base directory");
    let mut outfile = BufWriter::new(
        fs::File::create(THIS_BASEDIR.join(format!("tc_{filename}.out")))
            .expect("Could not create output file"),
    );
    for (mu, tc) in config.chempots().iter().zip(&tcs) {
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
    plot.set_domain(config.chempots().clone());
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
    plot.insert_image(tcsd2.iter().map(|(_, d2t)| *d2t).collect());
    plot.set_path(&format!(
        "{}/tcd2_{filename}.png",
        THIS_BASEDIR.to_string_lossy()
    ));
    plot.savefig().expect("Could not save figure");

    config
        .chempots()
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
    let (mq1, mq2, _mq3) = (0.350, 0.450, 1.5);
    let (mq1c, mq2c, _mq3c) = (0.125, 0.225, 1.2);

    // Adimensional
    let f0 = -0.876;

    // In units of m
    let fewtemps = (0., 0.15, 0.025);
    let moretemps = (0., 0.15, 0.01);
    let manytemps = (0., 0.15, 0.0001);

    // In GeV
    let fewchempots = (0., 0.8, 0.1);
    let morechempots = (0., 0.8, 0.025);
    let manychempots = (0., 0.8, 0.01);

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
        "nf = 2+1",
        "nf_2+1",
        "$n_{{F}}=2+1$",
    );

    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 ***");
    compute_propagators(&config);
    compute_transverse_propagators(&config);

    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1, MU AROUND THE QUARK MASSES ***");
    let oldfilename = config.filename; // save filename for later
    config.reset_temperatures((0., 0.05, 0.01)); // restrict to temperatures of interest
    config.reset_chemicalpotentials((mq1 - 0.05, mq1 + 0.05, 0.01));
    config.filename = "nf_2+1_m1"; // new filename
    compute_propagators(&config);
    config.reset_temperatures((0., 0.1, 0.02)); // restrict to temperatures of interest
    config.reset_chemicalpotentials((mq2 - 0.05, mq2 + 0.05, 0.01));
    config.filename = "nf_2+1_m2"; // new filename
    compute_propagators(&config);
    config.filename = oldfilename; // reset filename

    eprintln!("*** COMPUTING PROPAGATORS' IR LIMIT AT NF = 2 + 1 ***");
    config.reset_temperatures(moretemps);
    config.reset_chemicalpotentials(morechempots);
    compute_ir_limit(&config);

    eprintln!("*** COMPUTING PHASE DIAGRAM AT NF = 2 + 1 ***");
    config.reset_temperatures(manytemps);
    config.reset_chemicalpotentials(manychempots);
    let pb = compute_phase_diagram(&config);

    // Correct fieldconfig, set phase boundary and change filename for
    // calculations that follow
    config.reset_correctedfieldconfig(Some(&correctedfieldconfig));
    config.reset_phase_boundary(Some(&pb));
    config.filename = "nf_2+1_corrected";

    eprintln!("*** COMPUTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 ***");
    compute_phase_diagram(&config);

    eprintln!("*** FITTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 ***");
    let phase_boundary_params = parametrize_phase_boundary(
        config.phase_boundary().unwrap(),
        8,
        THIS_BASEDIR
            .join(format!("tc_{}_fit.png", config.filename))
            .to_str(),
    );
    let mut outfile = BufWriter::new(
        fs::File::create(THIS_BASEDIR.join(format!("tc_{}_fit_parameters.out", config.filename)))
            .expect("Could not create parameter file"),
    );
    writeln!(outfile, "{phase_boundary_params:?}").expect("Could not write to parameter file");

    eprintln!("*** COMPUTING CORRECTED PROPAGATORS AT NF = 2 + 1 ***");
    config.reset_temperatures(fewtemps);
    config.reset_chemicalpotentials(fewchempots);
    compute_propagators(&config);
    compute_transverse_propagators(&config);

    /* /* NF = 2 + 1 + 1 */
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
        "nf = 2+1+1",
        "nf_2+1+1",
        "$n_{{F}}=2+1+1$",
    );

    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 + 1 ***");
    compute_propagators(&config);

    eprintln!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 + 1, MU AROUND THE QUARK MASSES ***");
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

    eprintln!("*** COMPUTING PROPAGATORS' IR LIMIT AT NF = 2 + 1 + 1 ***");
    config.reset_temperatures(moretemps);
    config.reset_chemicalpotentials(morechempots);
    compute_ir_limit(&config);

    eprintln!("*** COMPUTING PHASE DIAGRAM AT NF = 2 + 1 + 1 ***");
    config.reset_temperatures(manytemps);
    // We already know the phase diagram thanks to the nf = 2 + 1 case,
    // let's go to larger chemical potentials to see what happens
    config.reset_chemicalpotentials((0., 3., 0.03));
    let pb = compute_phase_diagram(&config);

    config.reset_correctedfieldconfig(Some(&correctedfieldconfig));
    config.reset_phase_boundary(Some(&pb));
    config.filename = "nf_2+1+1_corrected";

    eprintln!("*** COMPUTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 + 1 ***");
    compute_phase_diagram(&config);

    eprintln!("*** COMPUTING CORRECTED PROPAGATORS AT NF = 2 + 1 + 1, MU AROUND THE LARGER QUARK MASS ***");
    let oldfilename = config.filename;
    config.reset_temperatures((0., 0.2, 0.025));
    config.reset_chemicalpotentials((mq3c - 0.2, mq3c + 1.3, 0.15));
    config.filename = "nf_2+1+1_corrected_m3";
    compute_propagators(&config);
    config.filename = oldfilename; */
}
