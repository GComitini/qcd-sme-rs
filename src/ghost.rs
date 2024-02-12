// The ghost functions are the same as those of pure Yang-Mills theory to one loop,
// so just re-export them. Unfortunately, cbindgen won't pick up the ffi functions,
// so we need to copy-paste those. While we are at it, we'll make them use this
// module's inlines rather than those of crate::ym.
//
// Do not write tests, it's not worth it at this stage. Keep in mind though that the
// ffi module is currently unchecked.

pub use crate::ym::ghost::native::*;

// For use in other languages, e.g. C/C++/Python
//
// - Re-export at crate::ffi, since symbols need to be unmangled anyway and the
//   namespace will not be preserved.
//
// - Here too each function calls its inline analogue, but the objective is to
//   switch from generic to concrete argument types so that the functions can
//   be compiled into a C dynamic library. To do so, we need to double their
//   number (one function for real arguments, another for complex arguments).
pub(crate) mod ffi {
    use super::inlines;
    use crate::{C, R};

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv_landau_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        ln_s_pl_1_2: R,
        g0: R,
    ) -> R {
        inlines::dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv_landau(s: R, g0: R) -> R {
        inlines::dressing_inv_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv_landau_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        ln_s_pl_1_2: C,
        g0: R,
    ) -> C {
        inlines::dressing_inv_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv_landau__complex(s: C, g0: R) -> C {
        inlines::dressing_inv_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        ln_s_pl_1_2: R,
        g0: R,
        xi: R,
    ) -> R {
        inlines::dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv(s: R, g0: R, xi: R) -> R {
        inlines::dressing_inv(s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        ln_s_pl_1_2: C,
        g0: R,
        xi: R,
    ) -> C {
        inlines::dressing_inv_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_inv__complex(s: C, g0: R, xi: R) -> C {
        inlines::dressing_inv(s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_landau_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        ln_s_pl_1_2: R,
        g0: R,
    ) -> R {
        inlines::dressing_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_landau(s: R, g0: R) -> R {
        inlines::dressing_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_landau_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        ln_s_pl_1_2: C,
        g0: R,
    ) -> C {
        inlines::dressing_landau_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_landau__complex(s: C, g0: R) -> C {
        inlines::dressing_landau(s, g0)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_sep(
        s: R,
        sinv: R,
        sinv2: R,
        s_pl_1_2: R,
        ln_s: R,
        ln_s_pl_1_2: R,
        g0: R,
        xi: R,
    ) -> R {
        inlines::dressing_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing(s: R, g0: R, xi: R) -> R {
        inlines::dressing(s, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing_sep__complex(
        s: C,
        sinv: C,
        sinv2: C,
        s_pl_1_2: C,
        ln_s: C,
        ln_s_pl_1_2: C,
        g0: R,
        xi: R,
    ) -> C {
        inlines::dressing_sep(s, sinv, sinv2, s_pl_1_2, ln_s, ln_s_pl_1_2, g0, xi)
    }

    #[no_mangle]
    pub extern "C" fn ghost__dressing__complex(s: C, g0: R, xi: R) -> C {
        inlines::dressing(s, g0, xi)
    }
}

pub(crate) mod inlines {
    #[allow(unused_imports)]
    pub use crate::ym::ghost::inlines::*;
}
