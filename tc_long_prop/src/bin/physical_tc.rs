use config::Config;
use lazy_static::lazy_static;
use log::info;
use serde::{Deserialize, Serialize};
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
        #[allow(clippy::too_many_arguments)]
        pub fn new(
            // (min, max, delta)
            temps: (R, R, R),
            chempots: (R, R, R),
            momenta: (R, R, R),
            om: R,
            renpoint: R,
            f0_init: R,
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

            // Needed because a different renormalization convention is used in the library
            let (nc, mg) = (fieldconfig.nc, fieldconfig.gluon);
            let f0 = f0_init
                + fieldconfig
                    .quarks
                    .iter()
                    .map(|(nf, mq)| (*nf as R) * (mg / mq).ln())
                    .sum::<R>()
                    * 4.
                    / (9. * (nc as R));

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

        pub fn set_temperatures(&mut self, temps: (R, R, R)) {
            self.temps = Self::compute_vec(temps);
        }

        pub fn set_chemicalpotentials(&mut self, chempots: (R, R, R)) {
            self.chempots = Self::compute_vec(chempots);
        }

        pub fn set_correctedfieldconfig(&mut self, correctedfieldconfig: Option<&'a FieldConfig>) {
            self.correctedfieldconfig = correctedfieldconfig;
            self.set_correctedf0();
        }

        pub fn set_phase_boundary(&mut self, phase_boundary: Option<&'a [(R, R)]>) {
            // If any, reduce the phase boundary before storing it
            self.phase_boundary = phase_boundary.map(reduce_phase_boundary);
        }

        // The boundary is taken to be part of the confined phase (inherited from
        // is_deconfined_phase)
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
                let cf0 = self.f0
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
                        / (9. * (self.fieldconfig.nc as R));
                info!("Set corrected f0 to {}", cf0);
                cf0
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

    fs::create_dir_all(THIS_BASEDIR.join("data")).expect("Could not create directory structure");

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
                info!("Computed (T/m, mu, p/m) = (0.0000, {mu:.4}, {p:.4}) for {label}. z = {z}.");
                val
            })
            .collect();
        let outfilename = format!(
            "{}_glulongprop_t_0.0000_mu_{mu:.4}.out",
            THIS_BASEDIR.join("data").join(filename).to_string_lossy()
        );
        let mut outfile = BufWriter::new(
            fs::File::create(&outfilename)
                .unwrap_or_else(|_| panic!("could not create {outfilename}")),
        );
        // Write pairs (p/m, m^2 Delta_L) to file
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
                // beta is dimensionful
                let beta = 1. / (t * m);
                let (fieldconfig, f0) = config.maybe_corrected_data(mu, *t);
                let z = propagator_l_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renfac;
                let outfilename = format!(
                    "{}_glulongprop_t_{:.4}_mu_{mu:.4}.out",
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
                        info!(
                            "Computed (T/m, mu, p/m) = ({t:.4}, {mu:.4}, {p:.4}) for {label}. z = {z}."
                        );
                        val
                    })
                    .collect();
                // Write pairs (p/m, m^2 Delta_L) to file
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

    fs::create_dir_all(THIS_BASEDIR.join("data")).expect("Could not create directory structure");

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
                info!(
                    "Computed (T/m, mu, p/m) = (0.0000, {mu:.4}, {p:.4}) for {label} (transverse). z = {z}."
                );
                val
            })
            .collect();
        let outfilename = format!(
            "{}_glutransprop_t_0.0000_mu_{mu:.4}.out",
            THIS_BASEDIR.join("data").join(filename).to_string_lossy()
        );
        let mut outfile = BufWriter::new(
            fs::File::create(&outfilename)
                .unwrap_or_else(|_| panic!("could not create {outfilename}")),
        );
        // Write pairs (p/m, m^2 Delta_T) to file
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
                // beta is dimensionful
                let beta = 1. / (t * m);
                let (fieldconfig, f0) = config.maybe_corrected_data(mu, *t);
                let z = propagator_t_landau(om, renpoint, beta, mu, f0, fieldconfig).re * renfac;
                let outfilename = format!(
                    "{}_glutransprop_t_{:.4}_mu_{mu:.4}.out",
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
                        let val = propagator_t_landau(om, p * m, beta, mu, f0, fieldconfig).re / z;
                        info!(
                            "Computed (T/m, mu, p/m) = ({t:.4}, {mu:.4}, {p:.4}) for {label} (transverse). z = {z}."
                        );
                        val
                    })
                    .collect();
                // Write pairs (p/m, m^2 Delta_T) to file
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
        info!("Computed (T/m, mu) = (0.0000, {mu:.4}) for {label}. z = {z}.");

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
                    info!("Computed (T/m, mu) = ({t:.4}, {mu:.4}) for {label}. z = {z}.");
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
        // Step 1: at fixed chemical potential, compute the propagator at p = p_min as a function of the temperature
        let (fieldconfig, f0) = config.maybe_corrected_data(mu, 0.);
        let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, fieldconfig).re * renfac;
        let mut vals =
            vec![propagator_l_zero_temp_landau(om, pbase * m, mu, f0, fieldconfig).re / z];
        info!("Computed (T/m, mu) = (0.0000, {mu:.4}) for {label}. z = {z}.");

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
                    info!("Computed (T/m, mu) = ({t:.4}, {mu:.4}) for {label}. z = {z}.");
                    val
                })
                .collect(),
        );
        // Step 2: at fixed chemical potential, find the critical temperature by maximizing the propagator
        let tc = find_critical_temperature(config.temps(), &vals);
        tcs.push(tc);
    }

    // Compute second derivative of critical temperature
    let (cpl, d2w) = (config.chempots().len(), D2W);
    let dmu = config.chempots()[1] - config.chempots()[0];
    let mut tcsd2 = Vec::with_capacity(cpl - 2 * d2w);
    let dmu2r = ((d2w * d2w) as R) * dmu * dmu;
    for i in d2w..cpl - d2w {
        let d2 = (tcs[i + d2w] + tcs[i - d2w] - 2. * tcs[i]) / dmu2r;
        tcsd2.push((config.chempots()[i], d2));
    }

    // Save critical temperature and second derivative to file
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

    // Plot critical temperature and second derivative
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

    // Return phase boundary
    config
        .chempots()
        .iter()
        .zip(tcs.iter())
        .map(|(&mu, &tc)| (mu, tc))
        .collect()
}

#[derive(Deserialize, Serialize)]
struct InputData {
    mg: f64,
    om: f64,
    renpoint: f64,
    mq: (f64, f64, f64),
    mqc: (f64, f64, f64),
    f0: f64,
    fewtemps: (f64, f64, f64),
    moretemps: (f64, f64, f64),
    manytemps: (f64, f64, f64),
    fewchempots: (f64, f64, f64),
    morechempots: (f64, f64, f64),
    manychempots: (f64, f64, f64),
    momenta: (f64, f64, f64),
}

fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    /* CRITICAL TEMPERATURE AS A FUNCTION OF THE CHEMICAL POTENTIAL
    FOR A PLAUSIBLE PHYSICAL CONFIGURATION, FIRST WITH NF = 2 + 1,
    THEN WITH NF = 2 + 1 + 1 */

    // In GeV
    let mut m = 0.656;
    let mut om = 1E-5;
    let mut renpoint = 4.;
    let (mut mq1, mut mq2, mut mq3) = (0.350, 0.450, 1.5);
    let (mut mq1c, mut mq2c, mut mq3c) = (0.125, 0.225, 1.2);

    // Adimensional
    let mut f0 = -0.876;

    // In units of m
    let mut fewtemps = (0., 0.1758, 0.0293);
    let mut moretemps = (0., 0.1758, 0.0117);
    let mut manytemps = (0., 0.1758, 0.00012);

    // In GeV
    let mut fewchempots = (0., 0.8, 0.1);
    let mut morechempots = (0., 0.8, 0.025);
    let mut manychempots = (0., 0.8, 0.01);

    // In units of m
    let mut momenta = (0.01, 3., 0.01);

    // Read and store configuration
    let configfilename = std::env::args().nth(1).unwrap_or("config.yml".to_string());
    if let Ok(configfile) = fs::File::open(&configfilename) {
        let inputdata: InputData =
            serde_yml::from_reader(configfile).expect("could not deserialize input data");
        info!(
            "Using configuration file {}",
            std::fs::canonicalize(&configfilename)
                .unwrap()
                .to_str()
                .unwrap()
        );
        (m, om, renpoint) = (inputdata.mg, inputdata.om, inputdata.renpoint);
        (mq1, mq2, mq3) = inputdata.mq;
        (mq1c, mq2c, mq3c) = inputdata.mqc;
        f0 = inputdata.f0;
        (fewtemps, moretemps, manytemps) =
            (inputdata.fewtemps, inputdata.moretemps, inputdata.manytemps);
        (fewchempots, morechempots, manychempots) = (
            inputdata.fewchempots,
            inputdata.morechempots,
            inputdata.manychempots,
        );
        momenta = inputdata.momenta;
    }

    let inputdata = InputData {
        mg: m,
        om,
        renpoint,
        mq: (mq1, mq2, mq3),
        mqc: (mq1c, mq2c, mq3c),
        f0,
        fewtemps,
        moretemps,
        manytemps,
        fewchempots,
        morechempots,
        manychempots,
        momenta,
    };

    fs::create_dir_all(THIS_BASEDIR.as_path()).expect("could not create base directory");
    let outconfigfile = BufWriter::new(
        fs::File::create(THIS_BASEDIR.join("config.yml"))
            .expect("could not create configuration file"),
    );
    serde_yml::to_writer(outconfigfile, &inputdata).expect("could not write to configuration file");

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

    info!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 ***");
    compute_propagators(&config);
    compute_transverse_propagators(&config);

    // * WARNING * This snippet needs to be re-written to account for non-hard-coded input data
    /* info!("*** COMPUTING PROPAGATORS AT NF = 2 + 1, MU AROUND THE QUARK MASSES ***");
    let oldfilename = config.filename; // save filename for later
    config.set_temperatures((0., 0.05, 0.01)); // restrict to temperatures of interest
    config.set_chemicalpotentials((mq1 - 0.05, mq1 + 0.05, 0.01));
    config.filename = "nf_2+1_m1"; // new filename
    compute_propagators(&config);
    config.set_temperatures((0., 0.1, 0.02)); // restrict to temperatures of interest
    config.set_chemicalpotentials((mq2 - 0.05, mq2 + 0.05, 0.01));
    config.filename = "nf_2+1_m2"; // new filename
    compute_propagators(&config);
    config.filename = oldfilename; // reset filename */

    info!("*** COMPUTING PROPAGATORS' IR LIMIT AT NF = 2 + 1 ***");
    config.set_temperatures(moretemps);
    config.set_chemicalpotentials(morechempots);
    compute_ir_limit(&config);

    info!("*** COMPUTING PHASE DIAGRAM AT NF = 2 + 1 ***");
    config.set_temperatures(manytemps);
    config.set_chemicalpotentials(manychempots);
    let pb = compute_phase_diagram(&config);

    // Correct fieldconfig, set phase boundary and change filename for
    // calculations that follow
    config.set_correctedfieldconfig(Some(&correctedfieldconfig));
    config.set_phase_boundary(Some(&pb));
    config.filename = "nf_2+1_corrected";

    info!("*** COMPUTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 ***");
    compute_phase_diagram(&config);

    info!("*** FITTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 ***");
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

    info!("*** COMPUTING CORRECTED PROPAGATORS AT NF = 2 + 1 ***");
    config.set_temperatures(fewtemps);
    config.set_chemicalpotentials(fewchempots);
    compute_propagators(&config);
    compute_transverse_propagators(&config);

    // * WARNING * What follows needs to be re-written to account for non-hard-coded input data
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

    info!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 + 1 ***");
    compute_propagators(&config);

    info!("*** COMPUTING PROPAGATORS AT NF = 2 + 1 + 1, MU AROUND THE QUARK MASSES ***");
    let oldfilename = config.filename;
    config.set_temperatures((0., 0.05, 0.01));
    config.set_chemicalpotentials((mq1 - 0.05, mq1 + 0.05, 0.01));
    config.filename = "nf_2+1+1_m1";
    compute_propagators(&config);
    config.set_temperatures((0., 0.1, 0.02));
    config.set_chemicalpotentials((mq2 - 0.05, mq2 + 0.05, 0.01));
    config.filename = "nf_2+1+1_m2";
    compute_propagators(&config);
    config.set_temperatures((0., 0.2, 0.04));
    config.set_chemicalpotentials((mq3 - 0.3, mq3 + 1.2, 0.15));
    config.filename = "nf_2+1+1_m3";
    compute_propagators(&config);
    config.filename = oldfilename;

    info!("*** COMPUTING PROPAGATORS' IR LIMIT AT NF = 2 + 1 + 1 ***");
    config.set_temperatures(moretemps);
    config.set_chemicalpotentials(morechempots);
    compute_ir_limit(&config);

    info!("*** COMPUTING PHASE DIAGRAM AT NF = 2 + 1 + 1 ***");
    config.set_temperatures(manytemps);
    // We already know the phase diagram thanks to the nf = 2 + 1 case,
    // let's go to larger chemical potentials to see what happens
    config.set_chemicalpotentials((0., 3., 0.03));
    let pb = compute_phase_diagram(&config);

    config.set_correctedfieldconfig(Some(&correctedfieldconfig));
    config.set_phase_boundary(Some(&pb));
    config.filename = "nf_2+1+1_corrected";

    info!("*** COMPUTING CORRECTED PHASE DIAGRAM AT NF = 2 + 1 + 1 ***");
    compute_phase_diagram(&config);

    info!("*** COMPUTING CORRECTED PROPAGATORS AT NF = 2 + 1 + 1, MU AROUND THE LARGER QUARK MASS ***");
    let oldfilename = config.filename;
    config.set_temperatures((0., 0.2, 0.025));
    config.set_chemicalpotentials((mq3c - 0.2, mq3c + 1.3, 0.15));
    config.filename = "nf_2+1+1_corrected_m3";
    compute_propagators(&config);
    config.filename = oldfilename; */
}
