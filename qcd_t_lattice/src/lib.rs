pub use qcd_sme::{qcd::FieldConfig, types::NCTYPE, C, R};
use std::path::Path;

pub const BASEDIR: &str = "target/qcd_t_lattice";

pub const OMEPS: R = 1E-6;
pub const NC: NCTYPE = 3;
pub const PBASE: R = 2E-3;
pub const PREN: R = 4.;

pub mod fit;

pub mod gluon {
    use qcd_sme::qcd::thermal::gluon::{
        propagator_l_landau_w_field_config, propagator_t_landau_w_field_config,
    };
    use qcd_sme::qcd::FieldConfig;
    use qcd_sme::{Num, C, R};

    pub fn propagator_l<T: Num>(om: T, p: R, beta: R, f0: R, config: &FieldConfig) -> C {
        propagator_l_landau_w_field_config(om, p, beta, 0., f0, config)
    }

    pub fn propagator_t<T: Num>(om: T, p: R, beta: R, f0: R, config: &FieldConfig) -> C {
        propagator_t_landau_w_field_config(om, p, beta, 0., f0, config)
    }
}

// Initialize env_logger
pub fn init<P: AsRef<Path>>(basedir: Option<P>) {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    if let Some(basedir) = &basedir {
        std::fs::create_dir_all(basedir)
            .unwrap_or_else(|e| panic!("Could not create base directory: {e}"));
    }
}
