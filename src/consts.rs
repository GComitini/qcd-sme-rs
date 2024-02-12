use crate::{types::NCTYPE, C};

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
