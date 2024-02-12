use crate::{
    types::{NCTYPE, NFTYPE},
    C, R,
};

/// The imaginary unit.
pub const I: C = C { re: 0., im: 1. };

/// The number of colors.
static mut NC: NCTYPE = 3;

/// Get the number of colors.
#[no_mangle]
pub extern "C" fn get_number_of_colors() -> NCTYPE {
    nc()
}

/// Set the number of colors.
#[no_mangle]
pub extern "C" fn set_number_of_colors(n: NCTYPE) {
    unsafe {
        NC = n;
    }
}

// This is just an alias to avoid typing "unsafe { NC }" every time.
#[inline(always)]
pub(crate) fn nc() -> NCTYPE {
    unsafe { NC }
}

/// The number of fermions.
static mut NF: NCTYPE = 1;

/// Get the number of fermions.
#[no_mangle]
pub extern "C" fn get_number_of_fermions() -> NFTYPE {
    nf()
}

/// Set the number of fermions.
#[no_mangle]
pub extern "C" fn set_number_of_fermions(n: NFTYPE) {
    unsafe {
        NF = n;
    }
}

// This is just an alias to avoid typing "unsafe { NF }" every time.
#[inline(always)]
pub(crate) fn nf() -> NCTYPE {
    unsafe { NF }
}

/// The default quark mass.
static mut M_QUARK: R = 0.3;

/// Get the default quark mass.
#[no_mangle]
pub extern "C" fn get_default_quark_mass() -> R {
    m_quark()
}

/// Set the number of colors.
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
