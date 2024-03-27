//! Constants, static variables and their getter/setter functions.
//!
//! Note: the static variables defined in this module are not multithread-safe.
//! Do not change them while multithreading!

use crate::{
    types::{NCTYPE, NFTYPE},
    C, R,
};

/// The imaginary unit.
pub const I: C = C { re: 0., im: 1. };

// The number of colors.
// If you change this line you MUST also change the value of NF_DIV_NC.
static mut NC: NCTYPE = 3;

/// Get the number of colors. Equals `3` if never changed.
#[no_mangle]
pub extern "C" fn get_number_of_colors() -> NCTYPE {
    nc()
}

/// Set the number of colors.
#[no_mangle]
pub extern "C" fn set_number_of_colors(n: NCTYPE) {
    unsafe {
        NC = n;
        NF_DIV_NC = (NF as R) / (NC as R);
    }
}

// This is just an alias to avoid typing "unsafe { NC }" every time.
#[inline(always)]
pub(crate) fn nc() -> NCTYPE {
    unsafe { NC }
}

// The number of fermions.
// If you change this line you MUST also change the value of NF_DIV_NC.
static mut NF: NCTYPE = 1;

/// Get the number of fermions. Equals `1` if never changed.
#[no_mangle]
pub extern "C" fn get_number_of_fermions() -> NFTYPE {
    nf()
}

/// Set the number of fermions.
#[no_mangle]
pub extern "C" fn set_number_of_fermions(n: NFTYPE) {
    unsafe {
        NF = n;
        NF_DIV_NC = (NF as R) / (NC as R);
    }
}

// This is just an alias to avoid typing "unsafe { NF }" every time.
#[inline(always)]
pub(crate) fn nf() -> NCTYPE {
    unsafe { NF }
}

/// The number of fermions divided by the number of colors.
static mut NF_DIV_NC: R = 1. / 3.;

// This is just an alias to avoid typing "unsafe { NF_DIV_NC }" every time.
#[inline(always)]
pub(crate) fn nf_div_nc() -> R {
    unsafe { NF_DIV_NC }
}

/// The default quark mass.
static mut M_QUARK: R = 0.5;

/// Get the default quark mass. Equals `0.5` (in arbitrary units) if never
/// changed.
#[no_mangle]
pub extern "C" fn get_default_quark_mass() -> R {
    m_quark()
}

/// Set the default quark mass.
#[no_mangle]
pub extern "C" fn set_default_quark_mass(m: R) {
    unsafe {
        M_QUARK = m;
    }
}

// This is just an alias to avoid typing "unsafe { M_QUARK }" every time.
#[inline(always)]
pub(crate) fn m_quark() -> R {
    unsafe { M_QUARK }
}

/// The default tolerance for numerical integrals.
static mut TOL_INTEGRAL: R = 1E-10;

/// Get the default tolerance for numerical integrals. Equals `1E-10` if never
/// changed.
#[no_mangle]
pub extern "C" fn get_default_tol_integral() -> R {
    tol_integral()
}

/// Set the default tolerance for numerical integrals.
///
/// # Panics
///
/// This function panics if `tol <= 0.`.
#[no_mangle]
pub extern "C" fn set_default_tol_integral(tol: R) {
    assert!(tol > 0.);
    unsafe {
        TOL_INTEGRAL = tol;
    }
}

// This is just an alias to avoid typing "unsafe { TOL_INTEGRAL }" every time.
#[inline(always)]
pub(crate) fn tol_integral() -> R {
    unsafe { TOL_INTEGRAL }
}

/// The default maximum number of iterations for numerical integrals.
static mut MAX_ITER_INTEGRAL: u32 = 50;

/// Get the default maximum number of iterations for numerical integrals.
/// Equals `50` if never changed.
#[no_mangle]
pub extern "C" fn get_default_max_iter_integral() -> u32 {
    max_iter_integral()
}

/// Set the default maximum number of iterations for numerical integrals.
#[no_mangle]
pub extern "C" fn set_default_max_iter_integral(iter: u32) {
    unsafe {
        MAX_ITER_INTEGRAL = iter;
    }
}

// This is just an alias to avoid typing "unsafe { MAX_ITER_INTEGRAL }" every time.
#[inline(always)]
pub(crate) fn max_iter_integral() -> u32 {
    unsafe { MAX_ITER_INTEGRAL }
}

/// Get the default integration method. Equals
/// `crate::Integral::G7K15(tol_integral(), max_iter_integral())` if never
/// changed.
pub fn get_default_integration_method() -> crate::Integral {
    crate::Integral::G7K15(tol_integral(), max_iter_integral())
}
