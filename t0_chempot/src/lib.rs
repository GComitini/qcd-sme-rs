use log::{debug, info};
use serde::{Deserialize, Serialize};
use std::io::{BufReader, BufWriter};
use std::path::Path;

pub use qcd_sme::{
    qcd::FieldConfig,
    types::{NCTYPE, NFTYPE},
    R,
};

use std::f64::consts::PI;

pub const BASEDIR: &str = "target/t0_chempot";

pub const OMEPS: R = 1E-6;

pub const NC_DEFAULT: NCTYPE = 3;
pub const F0_DEFAULT: R = -0.876;
pub const PBASE_DEFAULT: R = 2E-3;
pub const PREN_DEFAULT: R = 4.;
pub const MG_DEFAULT: R = 0.656;

#[derive(Debug, Copy, Clone)]
pub enum Component {
    T,
    L,
}

impl std::fmt::Display for Component {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl Component {
    pub fn abbrev(&self) -> &str {
        match self {
            Component::T => "trans.",
            Component::L => "long.",
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct InputConfig {
    pub fewchempots: [R; 3],
    pub manychempots: [R; 3],
    pub momenta: [R; 3],
    pub f0: R,
    pub pren: R,
    pub nc: NCTYPE,
    pub mg: R,
    pub quarkconfig: Vec<(NFTYPE, R)>,
    pub correctedquarkconfig: Option<Vec<(NFTYPE, R)>>,
    pub ymax: Option<R>,
}

impl Default for InputConfig {
    fn default() -> Self {
        Self {
            fewchempots: [0., 1., 0.1],
            manychempots: [0., 1., 0.001],
            momenta: [PBASE_DEFAULT, PREN_DEFAULT, 0.01],
            f0: F0_DEFAULT,
            pren: PREN_DEFAULT,
            nc: NC_DEFAULT,
            mg: MG_DEFAULT,
            quarkconfig: vec![(2, 0.35), (1, 0.45)],
            correctedquarkconfig: Some(vec![(2, 0.125), (1, 0.225)]),
            ymax: None,
        }
    }
}

impl InputConfig {
    pub fn load<P: AsRef<Path>>(path: P) -> Option<Self> {
        let Ok(file) = std::fs::File::open(&path) else {
            return None;
        };
        info!(
            "Using configuration file {}",
            std::fs::canonicalize(path).unwrap().to_str().unwrap()
        );
        let file = BufReader::new(file);
        serde_yml::from_reader(file)
            .unwrap_or_else(|e| panic!("Could not deserialize configuration file: {e}"))
    }
}

#[derive(Serialize, Debug, Clone)]
pub struct Config {
    #[serde(rename = "fewchempots")]
    fewchempots_init: [R; 3],
    #[serde(skip)]
    fewchempots: Vec<R>,
    #[serde(rename = "manychempots")]
    manychempots_init: [R; 3],
    #[serde(skip)]
    manychempots: Vec<R>,
    #[serde(rename = "momenta")]
    momenta_init: [R; 3],
    #[serde(skip)]
    momenta: Vec<R>,
    #[serde(rename = "f0")]
    f0_init: R,
    #[serde(skip)]
    f0: R,
    pub pren: R,
    pub nc: NCTYPE,
    pub mg: R,
    quarkconfig: Vec<(NFTYPE, R)>,
    correctedquarkconfig: Option<Vec<(NFTYPE, R)>>,
    #[serde(skip)]
    correctedf0: Option<R>,
    pub ymax: Option<R>,
}

impl Default for Config {
    fn default() -> Self {
        InputConfig::default().into()
    }
}

impl Config {
    pub fn load<P: AsRef<Path>>(p: P) -> Option<Self> {
        InputConfig::load(p).map(Into::into)
    }

    pub fn fewchempots_init(&self) -> [R; 3] {
        self.fewchempots_init
    }

    pub fn fewchempots(&self) -> &Vec<R> {
        &self.fewchempots
    }

    pub fn manychempots_init(&self) -> [R; 3] {
        self.manychempots_init
    }

    pub fn manychempots(&self) -> &Vec<R> {
        &self.manychempots
    }

    pub fn momenta(&self) -> &Vec<R> {
        &self.momenta
    }

    pub fn f0_init(&self) -> R {
        self.f0_init
    }

    pub fn fieldconfig(&self) -> FieldConfig {
        FieldConfig::new(self.nc, self.mg, self.quarkconfig.clone())
    }

    pub fn correctedfieldconfig(&self) -> Option<FieldConfig> {
        self.correctedquarkconfig
            .as_ref()
            .map(|cqc| FieldConfig::new(self.nc, self.mg, cqc.clone()))
    }

    pub fn min_quark_mass(&self) -> R {
        let mut mqm = R::INFINITY;
        for (_, mq) in self.quarkconfig.iter() {
            if *mq < mqm {
                mqm = *mq
            }
        }
        mqm
    }

    pub fn set_f0_init(&mut self, f0_init: R) {
        self.f0_init = f0_init;
        self.f0 = Self::compute_from_init_f0(f0_init, self.nc, &self.quarkconfig, self.mg);
        self.correctedf0 = Self::compute_corrected_f0(
            self.f0,
            self.nc,
            &self.quarkconfig,
            self.correctedquarkconfig.as_ref(),
        );
        if let Some(cf0) = self.correctedf0 {
            debug!("Set corrected f0 to {cf0}");
        }
    }

    pub fn maybe_corrected_data(&self, mu: R) -> (FieldConfig, R) {
        if mu <= self.min_quark_mass() || self.correctedquarkconfig.is_none() {
            debug!("Requested data correction. NOT correcting data for mu = {mu:.4}");
            (self.fieldconfig(), self.f0)
        } else {
            debug!("Requested data correction. CORRECTING data for mu = {mu:.4}");
            (
                self.correctedfieldconfig().unwrap(),
                self.correctedf0.unwrap(),
            )
        }
    }

    pub fn compute_from_init_f0(f0_init: R, nc: NCTYPE, qc: &[(NFTYPE, R)], mg: R) -> R {
        f0_init
            + qc.iter()
                .map(|(nf, mq)| (*nf as R) * (mg / mq).ln())
                .sum::<R>()
                * 4.
                / (9. * (nc as R))
    }

    pub fn compute_corrected_f0(
        f0: R,
        nc: NCTYPE,
        qc: &[(NFTYPE, R)],
        cqc: Option<&Vec<(NFTYPE, R)>>,
    ) -> Option<R> {
        cqc.map(|cqcu| {
            f0 + cqcu
                .iter()
                .enumerate()
                .map(|(i, &(nf, mqc))| {
                    let mq = qc[i].1;
                    (nf as R) * (mq / mqc).ln()
                })
                .sum::<R>()
                * 4.
                / (9. * (nc as R))
        })
    }

    fn compute_vec(value: [R; 3]) -> Vec<R> {
        let (min, max, delta) = (value[0], value[1], value[2]);
        let range = max - min;
        let n = (range / delta).round() as usize;
        (0..=n).map(|k| min + delta * (k as R)).collect()
    }
}

impl From<InputConfig> for Config {
    fn from(value: InputConfig) -> Self {
        let fewchempots: Vec<R> = Self::compute_vec(value.fewchempots)
            // Remove values of chemical potential around discontinuities
            .iter()
            .filter(|&&mu| {
                for qc in value.quarkconfig.iter() {
                    if mu > 0.99 * qc.1 && mu < 1.01 * qc.1 {
                        return false;
                    };
                }
                true
            })
            .copied()
            .collect();
        let manychempots = Self::compute_vec(value.manychempots)
            // Remove values of chemical potential around discontinuities
            .iter()
            .filter(|&&mu| {
                for qc in value.quarkconfig.iter() {
                    if mu > 0.99 * qc.1 && mu < 1.01 * qc.1 {
                        return false;
                    };
                }
                true
            })
            .copied()
            .collect();
        let momenta = Self::compute_vec(value.momenta);

        // Needed because a different renormalization convention is used in the library
        let (nc, mg) = (value.nc, value.mg);
        let f0_init = value.f0;
        let f0 = Self::compute_from_init_f0(f0_init, nc, &value.quarkconfig, mg);
        let cqc = &value.correctedquarkconfig;
        let correctedf0 = Self::compute_corrected_f0(f0, nc, &value.quarkconfig, cqc.as_ref());

        if let Some(cf0) = correctedf0 {
            debug!("Set corrected f0 to {cf0}");
        }

        Self {
            fewchempots_init: value.fewchempots,
            fewchempots,
            manychempots_init: value.manychempots,
            manychempots,
            momenta_init: value.momenta,
            momenta,
            f0_init,
            f0,
            pren: value.pren,
            nc,
            mg,
            quarkconfig: value.quarkconfig,
            correctedquarkconfig: value.correctedquarkconfig,
            correctedf0,
            ymax: value.ymax,
        }
    }
}

impl From<Config> for InputConfig {
    fn from(value: Config) -> Self {
        Self {
            fewchempots: value.fewchempots_init,
            manychempots: value.manychempots_init,
            momenta: value.momenta_init,
            f0: value.f0_init,
            pren: value.pren,
            nc: value.nc,
            mg: value.mg,
            quarkconfig: value.quarkconfig,
            correctedquarkconfig: value.correctedquarkconfig,
            ymax: value.ymax,
        }
    }
}

// Initialize env_logger, read configuration file and store it to result base directory
pub fn init<P: AsRef<Path>>(basedir: Option<P>) -> Config {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let config = Config::load("config.yml").unwrap_or_default();
    info!(
        "Configuration: nc = {}, mg = {}, f0_init = {}",
        config.nc, config.mg, config.f0_init
    );
    info!(
        "Configuration: quarkconfig = {:?}, corrected quarkconfig = {:?}",
        config.quarkconfig, config.correctedquarkconfig
    );
    if let Some(basedir) = &basedir {
        std::fs::create_dir_all(basedir)
            .unwrap_or_else(|e| panic!("Could not create base directory: {e}"));
        let outconfig = BufWriter::new(
            std::fs::File::create(basedir.as_ref().join("config.yml"))
                .unwrap_or_else(|e| panic!("Could not create configuration file: {e}")),
        );
        serde_yml::to_writer(outconfig, &config)
            .unwrap_or_else(|e| panic!("Could not deserialize configuration to file: {e}"));
    }
    config
}

pub mod gluon {
    use super::{Component, Config, R};
    use log::debug;

    pub use qcd_sme::qcd::thermal::gluon::{
        propagator_l_zero_temp_landau_w_field_config as propagator_l,
        propagator_t_zero_temp_landau_w_field_config as propagator_t,
    };

    pub fn propagator(om: R, p: R, mu: R, config: &Config, corr: bool, comp: Component) -> R {
        let pren = config.pren;
        let (fieldconfig, f0) = if corr {
            config.maybe_corrected_data(mu)
        } else {
            (config.fieldconfig(), config.f0)
        };
        let res = match comp {
            Component::L => {
                propagator_l(om, p, mu, f0, &fieldconfig).re
                    / (propagator_l(om, pren, mu, f0, &fieldconfig).re * pren * pren)
            }
            Component::T => {
                propagator_t(om, p, mu, f0, &fieldconfig).re
                    / (propagator_t(om, pren, mu, f0, &fieldconfig).re * pren * pren)
            }
        };
        debug!(
            "Computed{} {comp:?} propagator at (om, p, mu, f0) = ({om:.4}, {p:.4}, {mu:.4}, {f0:.4}): D_{comp:?} = {res:.5}", if corr {" (corrected)" } else {""}
        );
        res
    }

    pub(super) fn rawprop(p: R, config: &Config, mu: R, corr: bool, comp: Component) -> R {
        // Non-renormalized propagator at zero Matsubara frequency
        let (fieldconfig, f0) = if corr {
            config.maybe_corrected_data(mu)
        } else {
            (config.fieldconfig(), config.f0)
        };
        match comp {
            Component::L => propagator_l(super::OMEPS, p, mu, f0, &fieldconfig).re,
            Component::T => propagator_t(super::OMEPS, p, mu, f0, &fieldconfig).re,
        }
    }
}

pub fn alphastrong(config: &Config, mu: R, corr: bool, comp: Component) -> R {
    // This is alpha strong as implicitly defined by the parametrization we use
    // in the rest of the library. Although technically it depends on chemical
    // potential, on the component, etc., if the renormalization scale is large
    // wrt. the other scales (like the masses and mu) this dependence is very
    // small and can be neglected.
    // Note: this is NOT a running coupling, but a fixed-scale coupling.
    let pren = config.pren;
    4. * PI / (3. * (config.nc as R)) * pren * pren * gluon::rawprop(pren, config, mu, corr, comp)
}

pub fn f0strong(ast: R, config: &Config, mu: R, corr: bool, comp: Component) -> R {
    // This is f0_init in terms of alpha strong defined above
    let mut config = config.clone();
    let pren = config.pren;
    config.set_f0_init(0.);
    debug!("Note: the above \"Set corrected f0\" message refers to a dummy calculation");
    4. * PI / (3. * (config.nc as R) * ast)
        - 1. / (pren * pren * gluon::rawprop(pren, &config, mu, corr, comp))
}
