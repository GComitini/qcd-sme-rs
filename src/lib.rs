pub mod common;
pub mod consts;
pub mod ffi;
pub mod low_level;
pub mod types;
pub mod ym;

pub use consts::I;
pub use types::{Num, C, R};

/// The number of colors.
static mut NC: types::NCTYPE = 3;

/// Get the number of colors.
#[no_mangle]
pub extern "C" fn get_number_of_colors() -> types::NCTYPE {
    nc()
}

/// Set the number of colors.
#[no_mangle]
pub extern "C" fn set_number_of_colors(n: types::NCTYPE) {
    unsafe {
        NC = n;
    }
}

// This is just an alias to avoid typing "unsafe { NC }" every time
#[inline(always)]
pub(crate) fn nc() -> types::NCTYPE {
    unsafe { NC }
}
