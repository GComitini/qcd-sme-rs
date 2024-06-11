use peroxide::fuga::{LineStyle, Plot, Plot2D, PlotType, AD};
use qcd_sme::R;

pub mod prelude {
    pub use super::BASEDIR;
    pub use super::{find_critical_temperature, is_deconfined_phase, parametrize_phase_diagram};
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

pub trait ROrAD:
    std::ops::Add<Self, Output = Self> + std::ops::Mul<Self, Output = Self> + Sized + Copy
{
    const ZERO: Self;
}

impl ROrAD for R {
    const ZERO: R = 0.;
}

impl ROrAD for AD {
    const ZERO: AD = AD::AD0(0.);
}

fn parametric_phase_diagram<T: ROrAD>(mu: T, params: &[T]) -> T {
    let mut res = T::ZERO;
    for p in params[1..].iter().rev() {
        res = res + *p;
        res = res * mu;
    }
    res + params[0]
}

pub fn parametrize_phase_diagram(data: &[(R, R)], n: usize, plotpath: Option<&str>) -> Vec<R> {
    use peroxide::fuga::{matrix, Col, LevenbergMarquardt, Markers, Optimizer, AD1};

    let data_red = data.split(|(_, tc)| *tc == 0.).next().unwrap();
    let domain: Vec<R> = data_red.iter().map(|(mu, _)| *mu).collect();
    let values: Vec<R> = data_red.iter().map(|(_, tc)| *tc).collect();

    let params = Optimizer::new(
        peroxide::hstack!(domain.clone(), values.clone()),
        |mus, p| {
            Some(
                mus.iter()
                    .map(|mu| AD1(*mu, 0.))
                    .map(|mu| parametric_phase_diagram(mu, &p))
                    .collect(),
            )
        },
    )
    .set_init_param(vec![1.; n + 1])
    .set_max_iter(50)
    .set_method(LevenbergMarquardt)
    .set_lambda_init(1e-3)
    .set_lambda_max(1e+3)
    .optimize();

    if let Some(plotpath) = plotpath {
        let mut plot = Plot2D::new();
        plot.set_path(plotpath);
        plot.set_domain(domain.clone());
        plot.insert_image(values);
        plot.insert_image(
            domain
                .iter()
                .map(|mu| parametric_phase_diagram(*mu, &params))
                .collect(),
        );
        plot.set_plot_type(vec![(0, PlotType::Scatter)])
            .set_marker(vec![(0, Markers::Point)])
            .set_line_style(vec![(1, LineStyle::Solid)]);
        plot.savefig().expect("Could not save figure");
    }

    params
}

pub fn is_deconfined_phase(mu: R, t: R, phase_boundary: &[(R, R)]) -> bool {
    let phase_boundary = phase_boundary.split(|(_, t)| *t == 0.).next().unwrap();

    if t == 0. {
        return mu > phase_boundary.last().unwrap().0;
    }
    if mu == 0. {
        return t > phase_boundary.first().unwrap().1;
    }

    for (i, (mmu, tt)) in phase_boundary.iter().enumerate() {
        if *mmu == mu {
            return t > *tt;
        }
        if *tt == t {
            return mu > *mmu;
        }
        if *mmu > mu {
            // mmu is the smallest chemical potential in our set greater than mu,
            // therefore mu0 must be the largest chemical potential in our set less
            // than mu: mu0 < mu < mmu with no other chemical potentials in between.
            // Between these chemical potentials we interpolate the critical temperature
            // linearly.
            let (mu0, t0) = phase_boundary[i - 1];
            let tc = t0 + (tt - t0) / (mmu - mu0) * (mu - mu0);
            return t > tc;
        }
    }
    true
}
