use peroxide::fuga::{LineStyle, Plot, Plot2D, PlotType, AD};
use qcd_sme::R;

pub mod prelude {
    pub use super::BASEDIR;
    pub use super::{
        find_critical_temperature, is_deconfined_phase, parametrize_phase_boundary,
        reduce_phase_boundary,
    };
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
        propagator_t_landau_w_field_config as propagator_t_landau,
        propagator_t_zero_temp_landau_w_field_config as propagator_t_zero_temp_landau,
    };
    pub use qcd_sme::ym::gluon::propagator_landau as ym_propagator_l_zero_temp_landau;
    pub use qcd_sme::ym::thermal::gluon::propagator_l_landau as ym_propagator_l_landau;
}

pub const BASEDIR: &str = "target/tc_long_prop";

// Finds the global maximum (first argument) of its second argument
pub fn find_critical_temperature(ts: &[R], ds: &[R]) -> R {
    let (mut tc, mut dc) = (*ts.first().unwrap(), *ds.first().unwrap());
    for (&t, &d) in ts.iter().zip(ds.iter()).skip(1) {
        if d > dc {
            (tc, dc) = (t, d);
        }
    }
    tc
}

// Needed so we can define parametric_phase_diagram only once
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

// Sums the polynomial
fn parametric_phase_diagram<T: ROrAD>(mu: T, params: &[T]) -> T {
    let mut res = T::ZERO;
    for p in params[1..].iter().rev() {
        res = res + *p;
        res = res * mu;
    }
    res + params[0]
}

// This function selects the "physical" portion of the phase_boundary data (neglects
// data for mu > mu_c, after the critical temperature has vanished). The current
// implementation only works if the phase_boundary argument is fine enough to contain
// a chemical potential with Tc = 0. We assume we made it fine enough for this to be
// the case (i.e. we MUST make it fine enough). The implementation MUST be idempotent
// as this function is applied multiple times all throughout this crate
pub fn reduce_phase_boundary(phase_boundary: &[(R, R)]) -> &[(R, R)] {
    phase_boundary.split(|(_, t)| *t == 0.).next().unwrap()
}

// n is the degree of the polynomial used for fitting the phase boundary
pub fn parametrize_phase_boundary(pb: &[(R, R)], n: usize, plotpath: Option<&str>) -> Vec<R> {
    use peroxide::fuga::{matrix, Col, LevenbergMarquardt, Markers, Optimizer, AD1};

    // Data should already be reduced (i.e. points above mu = mu_c should have been removed)
    // at this point, but we do it again in case we want to reuse this function in code that
    // does not reduce them
    let pb_red = reduce_phase_boundary(pb);
    let domain: Vec<R> = pb_red.iter().map(|(mu, _)| *mu).collect();
    let values: Vec<R> = pb_red.iter().map(|(_, tc)| *tc).collect();

    // Fit the phase diagram to a polynomial of order n
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

// Determines whether the pair (mu, t) belongs to the deconfined phase given
// a phase boundary. The boundary itself is taken to be part of the confined
// phase for the sake of simplicity, and the algorithm assumes it has the
// expected shape. mu is assumed to be non-negative. Units must be passed
// with consistent dimensions and/or normalization.
pub fn is_deconfined_phase(mu: R, t: R, phase_boundary: &[(R, R)]) -> bool {
    // Again, we reduce the phase boundary even if it may not be needed
    let phase_boundary = reduce_phase_boundary(phase_boundary);

    // We do not extrapolate, the most we do is interpolate
    if mu < phase_boundary.first().unwrap().0 {
        panic!(
            "cannot extrapolate phase boundary data below the smallest provided chemical potential"
        );
    }

    if t == 0. {
        // This assumes the last chemical potential in phase_boundary is mu_c,
        // i.e. that the phase boundary is reduced and non-pathological (e.g.
        // it does not bend backwards)
        return mu > phase_boundary.last().unwrap().0;
    }

    if mu == 0. && phase_boundary.first().unwrap().0 == 0. {
        return t > phase_boundary.first().unwrap().1;
    }

    // We scan from lower to higher chemical potentials. If the chemical
    // potentials match we compare temperatures, if the temperatures match
    // we compare chemical potentials. Otherwise we interpolate (read below
    // and compare temperatures.
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
            // We interpolate the critical temperature between mu0 and mmu.
            let (mu0, t0) = phase_boundary[i - 1];
            let tc = t0 + (tt - t0) / (mmu - mu0) * (mu - mu0);
            return t > tc;
        }
    }
    true
}
