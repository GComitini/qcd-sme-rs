use qcd_sme::R;

pub mod prelude {
    pub use super::find_critical_temperature;
    pub use super::BASEDIR;
    pub use peroxide::util::plot::*;
    pub use qcd_sme::qcd::FieldConfig;
    pub use qcd_sme::{Num, C, R};
    pub use rayon::prelude::*;
}

pub mod one_mass {
    pub use super::prelude::*;
    pub use qcd_sme::consts::set_number_of_fermions;
    pub use qcd_sme::qcd::thermal::gluon::{propagator_l_landau, propagator_l_zero_temp_landau};
    pub use qcd_sme::ym::gluon::propagator_landau as ym_propagator_l_zero_temp_landau;
    pub use qcd_sme::ym::thermal::gluon::propagator_l_landau as ym_propagator_l_landau;
}

pub mod more_masses {
    pub use super::prelude::*;
    pub use qcd_sme::qcd::thermal::gluon::{
        propagator_l_landau_w_field_config as propagator_l_landau,
        propagator_l_zero_temp_landau_w_field_config as propagator_l_zero_temp_landau,
    };
    pub use qcd_sme::ym::gluon::propagator_landau as ym_propagator_l_zero_temp_landau;
    pub use qcd_sme::ym::thermal::gluon::propagator_l_landau as ym_propagator_l_landau;
}

pub const BASEDIR: &str = "target/tc_long_prop";

// Assumes the first pair to be computed at T=0
pub fn find_critical_temperature(ts: &[R], ds: &[R]) -> R {
    let first_d = *ds.first().unwrap();
    let (mut tc, mut dc) = (0., first_d);
    for (i, &d) in ds.iter().skip(1).enumerate() {
        if d > dc {
            tc = ts[i + 1];
            dc = d;
        }
    }
    tc
}
