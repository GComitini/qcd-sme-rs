pub mod ghost;
pub mod gluon;

pub use native::*;

// For use in Rust
//
// - Collect into a module to improve code organization, but immediately re-export.
//
// - Each function simply calls its inline analogue with the sole objective of not
//   inlining code when not necessary. The inline functions are still available to
//   the crate for more specialized use.
mod native {
    use super::inlines;
    use crate::{Num, C, R};

    pub fn tlog<T1: Num, T2: Num, T3: Copy + Into<C>>(q: R, om: T1, p: R, en: T2, delta: T3) -> C {
        inlines::tlog(q, om, p, en, delta)
    }

    pub fn tlog_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        inlines::tlog_same_mass(q, om, p, en)
    }

    pub fn tlog_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        inlines::tlog_zero_mass(q, om, p)
    }

    pub fn tlog_l<T1: Num, T2: Num, T3: Copy + Into<C>>(
        q: R,
        om: T1,
        p: R,
        en: T2,
        delta: T3,
    ) -> C {
        inlines::tlog_l(q, om, p, en, delta)
    }

    pub fn tlog_l_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        inlines::tlog_l_same_mass(q, om, p, en)
    }

    pub fn tlog_l_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        inlines::tlog_l_zero_mass(q, om, p)
    }

    pub fn tlog_t<T1: Num, T2: Num, T3: Copy + Into<C>>(
        q: R,
        om: T1,
        p: R,
        en: T2,
        delta: T3,
    ) -> C {
        inlines::tlog_t(q, om, p, en, delta)
    }

    pub fn tlog_t_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        inlines::tlog_t_same_mass(q, om, p, en)
    }

    pub fn tlog_t_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        inlines::tlog_t_zero_mass(q, om, p)
    }

    pub fn j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::j_m_i(q, m, beta)
    }

    pub fn j_0_i<T: Num>(q: R, beta: R) -> T {
        inlines::j_0_i(q, beta)
    }

    pub fn j_m_l_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::j_m_l_i(q, m, beta)
    }

    pub fn j_0_l_i<T: Num>(q: R, beta: R) -> T {
        inlines::j_0_l_i(q, beta)
    }

    pub fn j_m_t_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::j_m_t_i(q, m, beta)
    }

    pub fn j_0_t_i<T: Num>(q: R, beta: R) -> T {
        inlines::j_0_t_i(q, beta)
    }

    pub fn i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        inlines::i_m_m_i(q, om, p, m1, m2, beta)
    }

    pub fn i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::i_m_m_i_same_mass(q, om, p, m, beta)
    }

    pub fn i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::i_m_0_i(q, om, p, m, beta)
    }

    pub fn i_0_0_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        inlines::i_0_0_i(q, om, p, beta)
    }

    pub fn i_m_m_l_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        inlines::i_m_m_l_i(q, om, p, m1, m2, beta)
    }

    pub fn i_m_m_l_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::i_m_m_l_i_same_mass(q, om, p, m, beta)
    }

    pub fn i_m_0_l_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::i_m_0_l_i(q, om, p, m, beta)
    }

    pub fn i_0_0_l_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        inlines::i_0_0_l_i(q, om, p, beta)
    }

    pub fn i_m_m_t_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        inlines::i_m_m_t_i(q, om, p, m1, m2, beta)
    }

    pub fn i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::i_m_m_t_i_same_mass(q, om, p, m, beta)
    }

    pub fn i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::i_m_0_t_i(q, om, p, m, beta)
    }

    pub fn i_0_0_t_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        inlines::i_0_0_t_i(q, om, p, beta)
    }

    pub fn d_j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::d_j_m_i(q, m, beta)
    }

    pub fn d_j_0_i<T: Num>(q: R, beta: R) -> T {
        inlines::d_j_0_i(q, beta)
    }

    pub fn d2_j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::d2_j_m_i(q, m, beta)
    }

    pub fn d2_j_0_i<T: Num>(q: R, beta: R) -> T {
        inlines::d2_j_0_i(q, beta)
    }

    pub fn d_j_m_l_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::d_j_m_l_i(q, m, beta)
    }

    pub fn d_j_0_l_i<T: Num>(q: R, beta: R) -> T {
        inlines::d_j_0_l_i(q, beta)
    }

    pub fn d2_j_m_l_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::d2_j_m_l_i(q, m, beta)
    }

    pub fn d2_j_0_l_i<T: Num>(q: R, beta: R) -> T {
        inlines::d2_j_0_l_i(q, beta)
    }

    pub fn d_j_m_t_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::d_j_m_t_i(q, m, beta)
    }

    pub fn d_j_0_t_i<T: Num>(q: R, beta: R) -> T {
        inlines::d_j_0_t_i(q, beta)
    }

    pub fn d2_j_m_t_i<T: Num>(q: R, m: T, beta: R) -> T {
        inlines::d2_j_m_t_i(q, m, beta)
    }

    pub fn d2_j_0_t_i<T: Num>(q: R, beta: R) -> T {
        inlines::d2_j_0_t_i(q, beta)
    }

    pub fn d_i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        inlines::d_i_m_m_i(q, om, p, m1, m2, beta)
    }

    pub fn d_i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::d_i_m_m_i_same_mass(q, om, p, m, beta)
    }

    pub fn d_i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::d_i_m_0_i(q, om, p, m, beta)
    }

    pub fn d_i_0_0_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        inlines::d_i_0_0_i(q, om, p, beta)
    }

    pub fn d_i_m_m_t_i<T1: Num, T2: Num, T3: Num>(
        q: R,
        om: T1,
        p: R,
        m1: T2,
        m2: T3,
        beta: R,
    ) -> C {
        inlines::d_i_m_m_t_i(q, om, p, m1, m2, beta)
    }

    pub fn d_i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::d_i_m_m_t_i_same_mass(q, om, p, m, beta)
    }

    pub fn d_i_m_t_0_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, om, p, m, beta)
    }

    pub fn d_i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, om, p, m, beta)
    }

    pub fn d_i_0_0_t_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        inlines::d_i_0_0_t_i(q, om, p, beta)
    }
}

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
    pub extern "C" fn oneloop__tlog(q: R, om: R, p: R, en: R, delta: R) -> C {
        inlines::tlog(q, om, p, en, delta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_same_mass(q: R, om: R, p: R, en: R) -> C {
        inlines::tlog_same_mass(q, om, p, en)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_zero_mass(q: R, om: R, p: R) -> C {
        inlines::tlog_zero_mass(q, om, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_l(q: R, om: R, p: R, en: R, delta: R) -> C {
        inlines::tlog_l(q, om, p, en, delta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_l_same_mass(q: R, om: R, p: R, en: R) -> C {
        inlines::tlog_l_same_mass(q, om, p, en)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_l_zero_mass(q: R, om: R, p: R) -> C {
        inlines::tlog_l_zero_mass(q, om, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_t(q: R, om: R, p: R, en: R, delta: R) -> C {
        inlines::tlog_t(q, om, p, en, delta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_t_same_mass(q: R, om: R, p: R, en: R) -> C {
        inlines::tlog_t_same_mass(q, om, p, en)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__tlog_t_zero_mass(q: R, om: R, p: R) -> C {
        inlines::tlog_t_zero_mass(q, om, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__j_m_i(q: R, m: R, beta: R) -> R {
        inlines::j_m_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__j_0_i(q: R, beta: R) -> R {
        inlines::j_0_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__j_m_l_i(q: R, m: R, beta: R) -> R {
        inlines::j_m_l_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__j_0_l_i(q: R, beta: R) -> R {
        inlines::j_0_l_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__j_m_t_i(q: R, m: R, beta: R) -> R {
        inlines::j_m_t_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__j_0_t_i(q: R, beta: R) -> R {
        inlines::j_0_t_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_m_i(q: R, om: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_i(q, om, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_m_i_same_mass(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_m_i_same_mass(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_0_i(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_0_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_0_0_i(q: R, om: R, p: R, beta: R) -> C {
        inlines::i_0_0_i(q, om, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_m_l_i(q: R, om: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_l_i(q, om, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_m_l_i_same_mass(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_m_l_i_same_mass(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_0_l_i(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_0_l_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_0_0_l_i(q: R, om: R, p: R, beta: R) -> C {
        inlines::i_0_0_l_i(q, om, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_m_t_i(q: R, om: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_t_i(q, om, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_m_t_i_same_mass(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_m_t_i_same_mass(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_m_0_t_i(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_0_t_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__i_0_0_t_i(q: R, om: R, p: R, beta: R) -> C {
        inlines::i_0_0_t_i(q, om, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_j_m_i(q: R, m: R, beta: R) -> R {
        inlines::d_j_m_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_j_0_i(q: R, beta: R) -> R {
        inlines::d_j_0_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d2_j_m_i(q: R, m: R, beta: R) -> R {
        inlines::d2_j_m_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d2_j_0_i(q: R, beta: R) -> R {
        inlines::d2_j_0_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_j_m_l_i(q: R, m: R, beta: R) -> R {
        inlines::d_j_m_l_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_j_0_l_i(q: R, beta: R) -> R {
        inlines::d_j_0_l_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d2_j_m_l_i(q: R, m: R, beta: R) -> R {
        inlines::d2_j_m_l_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d2_j_0_l_i(q: R, beta: R) -> R {
        inlines::d2_j_0_l_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_j_m_t_i(q: R, m: R, beta: R) -> R {
        inlines::d_j_m_t_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_j_0_t_i(q: R, beta: R) -> R {
        inlines::d_j_0_t_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d2_j_m_t_i(q: R, m: R, beta: R) -> R {
        inlines::d2_j_m_t_i(q, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d2_j_0_t_i(q: R, beta: R) -> R {
        inlines::d2_j_0_t_i(q, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_m_m_i(q: R, om: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_i(q, om, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_m_m_i_same_mass(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::d_i_m_m_i_same_mass(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_m_0_i(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::d_i_m_0_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_0_0_i(q: R, om: R, p: R, beta: R) -> C {
        inlines::d_i_0_0_i(q, om, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_m_m_t_i(q: R, om: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_t_i(q, om, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_m_m_t_i_same_mass(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::d_i_m_m_t_i_same_mass(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_m_0_t_i(q: R, om: R, p: R, m: R, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, om, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__d_i_0_0_t_i(q: R, om: R, p: R, beta: R) -> C {
        inlines::d_i_0_0_t_i(q, om, p, beta)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use num::traits::Inv;

    use crate::common::{inlines::energy, thermal::inlines::bose_distribution_zero_chempot};
    use crate::{Num, C, I, R};
    use std::f64::consts::PI;

    const PI2: R = PI * PI;

    #[inline(always)]
    pub fn r_0<T1: Num, T2: Num, T3: Copy + Into<C>>(om: T1, p: R, en: T2, delta: T3) -> C {
        om * om + om * (en * 2. * I) + p * p + Into::<C>::into(delta)
    }

    #[inline(always)]
    pub fn r_0_same_mass<T1: Num, T2: Num>(om: T1, p: R, en: T2) -> C {
        om * om + om * (en * 2. * I) + p * p
    }

    #[inline(always)]
    pub fn tlog<T1: Num, T2: Num, T3: Copy + Into<C>>(q: R, om: T1, p: R, en: T2, delta: T3) -> C {
        let r_0 = r_0(om, p, en, delta);
        let q2p = q * 2. * p;
        ((r_0 + q2p) / (r_0 - q2p)).ln()
    }

    #[inline(always)]
    pub fn tlog_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        let r_0 = r_0_same_mass(om, p, en);
        let q2p = q * 2. * p;
        ((r_0 + q2p) / (r_0 - q2p)).ln()
    }

    #[inline(always)]
    pub fn tlog_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        let r_0 = r_0_same_mass(om, p, q);
        let q2p = q * 2. * p;
        ((r_0 + q2p) / (r_0 - q2p)).ln()
    }

    #[inline(always)]
    pub fn tlog_l<T1: Num, T2: Num, T3: Copy + Into<C>>(
        q: R,
        om: T1,
        p: R,
        en: T2,
        delta: T3,
    ) -> C {
        let ei2 = en * 2. * I;
        let pf = om * om + om * ei2 + Into::<C>::into(delta) + p * p * (om.inv() * ei2 + 1.);
        tlog(q, om, p, en, delta) * pf * pf
    }

    #[inline(always)]
    pub fn tlog_l_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        let ei2 = en * 2. * I;
        let pf = om * om + om * ei2 + p * p * (om.inv() * ei2 + 1.);
        tlog_same_mass(q, om, p, en) * pf * pf
    }

    #[inline(always)]
    pub fn tlog_l_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        let ei2 = q * 2. * I;
        let pf = om * om + om * ei2 + p * p * (om.inv() * ei2 + 1.);
        tlog_zero_mass(q, om, p) * pf * pf
    }

    #[inline(always)]
    pub fn tlog_t<T1: Num, T2: Num, T3: Copy + Into<C>>(
        q: R,
        om: T1,
        p: R,
        en: T2,
        delta: T3,
    ) -> C {
        let r0 = r_0(om, p, en, delta);
        let q2p = q * 2. * p;
        (r0 + q2p) * (r0 - q2p) * tlog(q, om, p, en, delta)
    }

    #[inline(always)]
    pub fn tlog_t_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        let r0 = r_0_same_mass(om, p, en);
        let q2p = q * 2. * p;
        (r0 + q2p) * (r0 - q2p) * tlog_same_mass(q, om, p, en)
    }

    #[inline(always)]
    pub fn tlog_t_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        let r0 = r_0_same_mass(om, p, q);
        let q2p = q * 2. * p;
        (r0 + q2p) * (r0 - q2p) * tlog_zero_mass(q, om, p)
    }

    #[inline(always)]
    pub fn j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        let q2 = q * q;
        let en = energy(q, m);
        bose_distribution_zero_chempot(en, beta) * q2 / (en * 2. * PI2)
    }

    #[inline(always)]
    pub fn j_0_i<T: Num>(q: R, beta: R) -> T {
        Into::<T>::into(bose_distribution_zero_chempot(q, beta) * q / (2. * PI2))
    }

    #[inline(always)]
    pub fn j_m_l_i<T: Num>(q: R, m: T, beta: R) -> T {
        let q2 = q * q;
        let en = energy(q, m);
        -bose_distribution_zero_chempot(en, beta) * q2 * en / (2. * PI2)
    }

    #[inline(always)]
    pub fn j_0_l_i<T: Num>(q: R, beta: R) -> T {
        let q2 = q * q;
        let q3 = q2 * q;
        Into::<T>::into(-bose_distribution_zero_chempot(q, beta) * q3 / (2. * PI2))
    }

    #[inline(always)]
    pub fn j_m_t_i<T: Num>(q: R, m: T, beta: R) -> T {
        let q2 = q * q;
        let q4 = q2 * q2;
        let en = energy(q, m);
        bose_distribution_zero_chempot(en, beta) * q4 / (en * 6. * PI2)
    }

    #[inline(always)]
    pub fn j_0_t_i<T: Num>(q: R, beta: R) -> T {
        let q2 = q * q;
        let q3 = q2 * q;
        Into::<T>::into(bose_distribution_zero_chempot(q, beta) * q3 / (6. * PI2))
    }

    #[inline(always)]
    pub fn i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1
            * (tlog(q, om, p, en1, delta) + tlog(q, -om, p, en1, delta));
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2
            * (tlog(q, om, p, en2, -delta) + tlog(q, -om, p, en2, -delta));
        (t1 + t2) * q / (p * 16. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let t = bose_distribution_zero_chempot(en, beta) / en
            * (tlog_same_mass(q, om, p, en) + tlog_same_mass(q, -om, p, en));
        t * q / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let delta = -m * m;
        let t1 = bose_distribution_zero_chempot(en, beta) * q / en
            * (tlog(q, om, p, en, delta) + tlog(q, -om, p, en, delta));
        let t2 = bose_distribution_zero_chempot(q, beta)
            * (tlog(q, om, p, q, -delta) + tlog(q, -om, p, q, -delta));
        (t1 + t2) / (p * 16. * PI2)
    }

    #[inline(always)]
    pub fn i_0_0_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        let t = bose_distribution_zero_chempot(q, beta)
            * (tlog_zero_mass(q, om, p) + tlog_zero_mass(q, -om, p));
        t / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_l_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1
            * (tlog_l(q, om, p, en1, delta) + tlog_l(q, -om, p, en1, delta)
                - 8. * qp * (om2_plus_p2 + delta));
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2
            * (tlog_l(q, om, p, en2, -delta) + tlog_l(q, -om, p, en2, -delta)
                - 8. * qp * (om2_plus_p2 - delta));
        (om2_plus_p2 * p3 * 64. * PI2).inv() * om2 * q * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_m_m_l_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t = bose_distribution_zero_chempot(en, beta) / en
            * (-om2_plus_p2 * 8. * qp
                + tlog_l_same_mass(q, om, p, en)
                + tlog_l_same_mass(q, -om, p, en));
        (om2_plus_p2 * p3 * 32. * PI2).inv() * om2 * q * t
    }

    #[inline(always)]
    pub fn i_m_0_l_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let delta = -Into::<C>::into(m * m);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t1 = bose_distribution_zero_chempot(en, beta) * q / en
            * (tlog_l(q, om, p, en, delta) + tlog_l(q, -om, p, en, delta)
                - 8. * qp * (om2_plus_p2 + delta));
        let t2 = bose_distribution_zero_chempot(q, beta)
            * (tlog_l(q, om, p, q, -delta) + tlog_l(q, -om, p, q, -delta)
                - 8. * qp * (om2_plus_p2 - delta));
        (om2_plus_p2 * p3 * 64. * PI2).inv() * om2 * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_0_0_l_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t = bose_distribution_zero_chempot(q, beta)
            * (-om2_plus_p2 * 8. * qp + tlog_l_zero_mass(q, om, p) + tlog_l_zero_mass(q, -om, p));
        (om2_plus_p2 * p3 * 32. * PI2).inv() * om2 * t
    }

    #[inline(always)]
    pub fn i_m_m_t_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1
            * (tlog_t(q, om, p, en1, delta) + tlog_t(q, -om, p, en1, delta)
                - 8. * qp * (om2_plus_p2 + delta));
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2
            * (tlog_t(q, om, p, en2, -delta) + tlog_t(q, -om, p, en2, -delta)
                - 8. * qp * (om2_plus_p2 - delta));
        -(p3 * 128. * PI2).inv() * q * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t = bose_distribution_zero_chempot(en, beta) / en
            * (-om2_plus_p2 * 8. * qp
                + tlog_t_same_mass(q, om, p, en)
                + tlog_t_same_mass(q, -om, p, en));
        -(p3 * 64. * PI2).inv() * q * t
    }

    #[inline(always)]
    pub fn i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let delta = -Into::<C>::into(m * m);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t1 = bose_distribution_zero_chempot(en, beta) * q / en
            * (tlog_t(q, om, p, en, delta) + tlog_t(q, -om, p, en, delta)
                - 8. * qp * (om2_plus_p2 + delta));
        let t2 = bose_distribution_zero_chempot(q, beta)
            * (tlog_t(q, om, p, q, -delta) + tlog_t(q, -om, p, q, -delta)
                - 8. * qp * (om2_plus_p2 - delta));
        -(p3 * 128. * PI2).inv() * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_0_0_t_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let om2 = om * om;
        let om2_plus_p2 = om2 + p2;
        let t = bose_distribution_zero_chempot(q, beta)
            * (-om2_plus_p2 * 8. * qp + tlog_t_zero_mass(q, om, p) + tlog_t_zero_mass(q, -om, p));
        -(p3 * 64. * PI2).inv() * t
    }

    #[inline(always)]
    pub fn d_j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        -bose_distribution_zero_chempot(en, beta) / (en * 4. * PI2)
    }

    #[inline(always)]
    pub fn d_j_0_i<T: Num>(q: R, beta: R) -> T {
        Into::<T>::into(-bose_distribution_zero_chempot(q, beta) / (q * 4. * PI2))
    }

    #[inline(always)]
    pub fn d2_j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        let en3 = en * en * en;
        let ben = en * beta;
        let bp = bose_distribution_zero_chempot(en, beta);
        let bm = bose_distribution_zero_chempot(-en, beta);
        bp / (en3 * 8. * PI2) * (-ben * bm + 1.)
    }

    #[inline(always)]
    pub fn d2_j_0_i<T: Num>(q: R, beta: R) -> T {
        let q3 = q * q * q;
        let bq = q * beta;
        let bp = bose_distribution_zero_chempot(q, beta);
        let bm = bose_distribution_zero_chempot(-q, beta);
        Into::<T>::into(bp / (q3 * 8. * PI2) * (-bq * bm + 1.))
    }

    #[inline(always)]
    pub fn d_j_m_l_i<T: Num>(q: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        en * bose_distribution_zero_chempot(en, beta) / (4. * PI2)
    }

    #[inline(always)]
    pub fn d_j_0_l_i<T: Num>(q: R, beta: R) -> T {
        Into::<T>::into(q * bose_distribution_zero_chempot(q, beta) / (4. * PI2))
    }

    #[inline(always)]
    pub fn d2_j_m_l_i<T: Num>(q: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        let bp = bose_distribution_zero_chempot(en, beta);
        let bm = bose_distribution_zero_chempot(-en, beta);
        bp / (8. * PI2) * (bm * beta + en.inv())
    }

    #[inline(always)]
    pub fn d2_j_0_l_i<T: Num>(q: R, beta: R) -> T {
        let bp = bose_distribution_zero_chempot(q, beta);
        let bm = bose_distribution_zero_chempot(-q, beta);
        Into::<T>::into(bp / (8. * PI2) * (bm * beta + 1. / q))
    }

    #[inline(always)]
    pub fn d_j_m_t_i<T: Num>(q: R, m: T, beta: R) -> T {
        let q2 = q * q;
        let en = energy(q, m);
        -bose_distribution_zero_chempot(en, beta) * q2 / (en * 4. * PI2)
    }

    #[inline(always)]
    pub fn d_j_0_t_i<T: Num>(q: R, beta: R) -> T {
        Into::<T>::into(-bose_distribution_zero_chempot(q, beta) * q / (4. * PI2))
    }

    #[inline(always)]
    pub fn d2_j_m_t_i<T: Num>(q: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        bose_distribution_zero_chempot(en, beta) / (en * 8. * PI2)
    }

    #[inline(always)]
    pub fn d2_j_0_t_i<T: Num>(q: R, beta: R) -> T {
        Into::<T>::into(bose_distribution_zero_chempot(q, beta) / (q * 8. * PI2))
    }

    #[inline(always)]
    pub fn d_i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, p: R, m1: T2, m2: T3, beta: R) -> C {
        let q2p = q * p * 2.;
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let r01 = r_0(om, p, en1, delta);
        let r02 = r_0(om, p, en2, -delta);
        let r1pinv = (r01 + q2p).inv();
        let r1minv = (r01 - q2p).inv();
        let r2pinv = (r02 + q2p).inv();
        let r2minv = (r02 - q2p).inv();
        let r01_opp = r_0(-om, p, en1, delta);
        let r02_opp = r_0(-om, p, en2, -delta);
        let r1pinv_opp = (r01_opp + q2p).inv();
        let r1minv_opp = (r01_opp - q2p).inv();
        let r2pinv_opp = (r02_opp + q2p).inv();
        let r2minv_opp = (r02_opp - q2p).inv();
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let t10 = b1 / en1 * (r1pinv + r1minv + r1pinv_opp + r1minv_opp);
        let t1 = b1 / en1 * (r1pinv - r1minv + r1pinv_opp - r1minv_opp);
        let t2 = b2 / en2 * (r2pinv - r2minv + r2pinv_opp - r2minv_opp);
        (-t10 + (t2 - t1) * q / p) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let q2p = q * p * 2.;
        let en = energy(q, m);
        let r0 = r_0_same_mass(om, p, en);
        let rpinv = (r0 + q2p).inv();
        let rminv = (r0 - q2p).inv();
        let r0_opp = r_0_same_mass(-om, p, en);
        let rpinv_opp = (r0_opp + q2p).inv();
        let rminv_opp = (r0_opp - q2p).inv();
        -bose_distribution_zero_chempot(en, beta)
            * (en * 8. * PI2).inv()
            * (rpinv + rminv + rpinv_opp + rminv_opp)
    }

    #[inline(always)]
    pub fn d_i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let q2p = q * p * 2.;
        let delta = -m * m;
        let en = energy(q, m);
        let r01 = r_0(om, p, en, delta);
        let r02 = r_0(om, p, q, -delta);
        let r1pinv = (r01 + q2p).inv();
        let r1minv = (r01 - q2p).inv();
        let r2pinv = (r02 + q2p).inv();
        let r2minv = (r02 - q2p).inv();
        let r01_opp = r_0(-om, p, en, delta);
        let r02_opp = r_0(-om, p, q, -delta);
        let r1pinv_opp = (r01_opp + q2p).inv();
        let r1minv_opp = (r01_opp - q2p).inv();
        let r2pinv_opp = (r02_opp + q2p).inv();
        let r2minv_opp = (r02_opp - q2p).inv();
        let b1 = bose_distribution_zero_chempot(en, beta);
        let b2 = bose_distribution_zero_chempot(q, beta);
        let t10 = b1 / en * (r1pinv + r1minv + r1pinv_opp + r1minv_opp);
        let t1 = b1 / en * (r1pinv - r1minv + r1pinv_opp - r1minv_opp);
        let t2 = b2 / q * (r2pinv - r2minv + r2pinv_opp - r2minv_opp);
        (-t10 + (t2 - t1) * q / p) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        let q2p = q * p * 2.;
        let r0 = r_0_same_mass(om, p, q);
        let rpinv = (r0 + q2p).inv();
        let rminv = (r0 - q2p).inv();
        let r0_opp = r_0_same_mass(-om, p, q);
        let rpinv_opp = (r0_opp + q2p).inv();
        let rminv_opp = (r0_opp - q2p).inv();
        -bose_distribution_zero_chempot(q, beta) * (rpinv + rminv + rpinv_opp + rminv_opp)
            / (q * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_t_i<T1: Num, T2: Num, T3: Num>(
        q: R,
        om: T1,
        p: R,
        m1: T2,
        m2: T3,
        beta: R,
    ) -> C {
        let p2 = p * p;
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let tlog1p = tlog(q, om, p, en1, delta);
        let tlog1m = tlog(q, -om, p, en1, delta);
        let r01p = r_0(om, p, en1, delta);
        let r01m = r_0(-om, p, en1, delta);
        let en2 = energy(q, m2);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let tlog2p = tlog(q, om, p, en2, -delta);
        let tlog2m = tlog(q, -om, p, en2, -delta);
        let r02p = r_0(om, p, en2, -delta);
        let r02m = r_0(-om, p, en2, -delta);
        let t1 = -b1 / en1 * (tlog1p + tlog1m) / 4.;
        let t2 = (b2 / en2 - Into::<C>::into(b1 / en1)) * q / p;
        let t3 = (b1 / en1 * (r01p * tlog1p + r01m * tlog1m)
            - b2 / en2 * (r02p * tlog2p + r02m * tlog2m))
            / (p2 * 8.);
        (t1 + t2 + t3) * q / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let b = bose_distribution_zero_chempot(en, beta);
        let tlogp = tlog_same_mass(q, om, p, en);
        let tlogm = tlog_same_mass(q, -om, p, en);
        -(en * p * 32. * PI2).inv() * b * (tlogp + tlogm) * q
    }

    #[inline(always)]
    pub fn d_i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, p: R, m: T2, beta: R) -> C {
        let p2 = p * p;
        let delta = -m * m;
        let en1 = energy(q, m);
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let tlog1p = tlog(q, om, p, en1, delta);
        let tlog1m = tlog(q, -om, p, en1, delta);
        let r01p = r_0(om, p, en1, delta);
        let r01m = r_0(-om, p, en1, delta);
        let b2 = bose_distribution_zero_chempot(q, beta);
        let tlog2p = tlog(q, om, p, q, -delta);
        let tlog2m = tlog(q, -om, p, q, -delta);
        let r02p = r_0(om, p, q, -delta);
        let r02m = r_0(-om, p, q, -delta);
        let t1 = -b1 / en1 * (tlog1p + tlog1m) / 4.;
        let t2 = (b2 / q - Into::<C>::into(b1 / en1)) * q / p;
        let t3 = (b1 / en1 * (r01p * tlog1p + r01m * tlog1m)
            - b2 / q * (r02p * tlog2p + r02m * tlog2m))
            / (p2 * 8.);
        (t1 + t2 + t3) * q / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_t_i<T: Num>(q: R, om: T, p: R, beta: R) -> C {
        -bose_distribution_zero_chempot(q, beta)
            * (tlog_same_mass(q, om, p, q) + tlog_same_mass(q, -om, p, q))
            / (p * 32. * PI2)
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::thermal;
    use crate::{Num, C, I, R};

    const TOLERANCE: R = 1e-11;

    fn assert_equal<T: Num>(lhs: T, rhs: T) {
        if rhs != T::zero() {
            assert!(
                (lhs / rhs - 1.).abs() < TOLERANCE,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e} rhs"
            );
        } else {
            assert!(
                (lhs - rhs).abs() < TOLERANCE,
                "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
            );
        }
    }

    #[test]
    fn test_tlog() {
        use thermal::{ffi::oneloop__tlog, tlog};

        let args: [(R, R, R, R, R); 6] = [
            (0.03, 0., 2.12, 4., 2.),
            (0.047, 0., 2.12, 4., 2.),
            (0.047, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 1.27, 4., 2.),
            (0.047, 0.18, 1.27, 6.9, 2.),
            (0.047, 0.18, 1.27, 6.9, 3.21),
        ];

        let res: [C; 6] = [
            0.039177220079546576 + 0. * I,
            0.06138906758007032 + 0. * I,
            0.05824470539037811 - 0.012858080275916684 * I,
            0.05666490383705108 - 0.022398146165630385 * I,
            0.0447257555087438 - 0.03049215271192196 * I,
            0.03897509874778243 - 0.01994626073948082 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(tlog(*q, *om, *p, *en, *delta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(oneloop__tlog(*q, *om, *p, *en, *delta), res[i])
            });
    }

    #[test]
    fn test_tlog_same_mass() {
        use thermal::{ffi::oneloop__tlog_same_mass, tlog_same_mass};

        let args: [(R, R, R, R); 5] = [
            (0.03, 0., 2.12, 4.),
            (0.047, 0., 2.12, 4.),
            (0.047, 0.18, 2.12, 4.),
            (0.047, 0.18, 1.27, 4.),
            (0.047, 0.18, 1.27, 6.9),
        ];

        let res: [C; 5] = [
            0.05661889399950811 + 0. * I,
            0.08873742845995349 + 0. * I,
            0.07998357415820304 - 0.025473106952291366 * I,
            0.0821116772242251 - 0.07200880785684338 * I,
            0.04420928411087925 - 0.06681663995758536 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en))| assert_equal(tlog_same_mass(*q, *om, *p, *en), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p, en))| {
            assert_equal(oneloop__tlog_same_mass(*q, *om, *p, *en), res[i])
        });
    }

    #[test]
    fn test_tlog_zero_mass() {
        use thermal::{ffi::oneloop__tlog_zero_mass, tlog_zero_mass};

        let args: [(R, R, R); 4] = [
            (0.03, 0., 2.12),
            (0.047, 0., 2.12),
            (0.047, 0.18, 2.12),
            (0.047, 0.18, 1.27),
        ];

        let res: [C; 4] = [
            0.05661889399950811 + 0. * I,
            0.08873742845995349 + 0. * I,
            0.08810024178455797 - 0.00032972192107351254 * I,
            0.1453563556448732 - 0.0015000913395812037 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p))| assert_equal(tlog_zero_mass(*q, *om, *p), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p))| assert_equal(oneloop__tlog_zero_mass(*q, *om, *p), res[i]));
    }

    #[test]
    fn test_tlog_l() {
        use thermal::{ffi::oneloop__tlog_l, tlog_l};

        let args: [(R, R, R, R, R); 6] = [
            (0.03, 0.001, 2.12, 4., 2.),
            (0.047, 0.001, 2.12, 4., 2.),
            (0.047, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 1.27, 4., 2.),
            (0.047, 0.18, 1.27, 6.9, 2.),
            (0.047, 0.18, 1.27, 6.9, 3.21),
        ];

        let res: [C; 6] = [
            -50647299.45180459 + 80701.18909865098 * I,
            -79362202.76459928 + 126491.79788072132 * I,
            -2321.371274775795 + 672.8860941248199 * I,
            -290.303789411238 + 149.6786403255582 * I,
            -683.0051933688675 + 525.8933457752748 * I,
            -594.7903868483935 + 364.63932196183583 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(tlog_l(*q, *om, *p, *en, *delta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(oneloop__tlog_l(*q, *om, *p, *en, *delta), res[i])
            });
    }

    #[test]
    fn test_tlog_l_same_mass() {
        use thermal::{ffi::oneloop__tlog_l_same_mass, tlog_l_same_mass};

        let args: [(R, R, R, R); 5] = [
            (0.03, 0.001, 2.12, 4.),
            (0.047, 0.001, 2.12, 4.),
            (0.047, 0.18, 2.12, 4.),
            (0.047, 0.18, 1.27, 4.),
            (0.047, 0.18, 1.27, 6.9),
        ];

        let res: [C; 5] = [
            -73195320.67949668 + 148655.6384399857 * I,
            -114717261.85563211 + 233143.37255073342 * I,
            -3189.5257931497435 + 1176.2653331677252 * I,
            -421.5170034292345 + 404.6074260300294 * I,
            -675.5698424879727 + 1081.3032612270954 * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, en))| {
            assert_equal(tlog_l_same_mass(*q, *om, *p, *en), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, en))| {
            assert_equal(oneloop__tlog_l_same_mass(*q, *om, *p, *en), res[i])
        });
    }

    #[test]
    fn test_tlog_l_zero_mass() {
        use thermal::{ffi::oneloop__tlog_l_zero_mass, tlog_l_zero_mass};

        let args: [(R, R, R); 4] = [
            (0.03, 0.001, 2.12),
            (0.047, 0.001, 2.12),
            (0.047, 0.18, 2.12),
            (0.047, 0.18, 1.27),
        ];

        let res: [C; 4] = [
            -4116.1062331546 + 137.29670408557539 * I,
            -15836.410087764021 + 337.31483528723743 * I,
            1.32005388356121 + 1.8806664241476745 * I,
            0.29041395973279666 + 0.4080161453514444 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p))| assert_equal(tlog_l_zero_mass(*q, *om, *p), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p))| {
            assert_equal(oneloop__tlog_l_zero_mass(*q, *om, *p), res[i])
        });
    }

    #[test]
    fn test_tlog_t() {
        use thermal::{ffi::oneloop__tlog_t, tlog_t};

        let args: [(R, R, R, R, R); 6] = [
            (0.03, 0.001, 2.12, 4., 2.),
            (0.047, 0.001, 2.12, 4., 2.),
            (0.047, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 1.27, 4., 2.),
            (0.047, 0.18, 1.27, 6.9, 2.),
            (0.047, 0.18, 1.27, 6.9, 3.21),
        ];

        let res: [C; 6] = [
            1.651753049290665 + 0.00203572060840294 * I,
            2.586783393101998 + 0.0031904825642123558 * I,
            2.599779527579854 + 0.5742667042452293 * I,
            0.8698134821767374 + 0.34402714128455986 * I,
            0.8699268749106882 + 0.593369498685203 * I,
            1.1588811274128774 + 0.5932693242723573 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(tlog_t(*q, *om, *p, *en, *delta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal(oneloop__tlog_t(*q, *om, *p, *en, *delta), res[i])
            });
    }

    #[test]
    fn test_tlog_t_same_mass() {
        use thermal::{ffi::oneloop__tlog_t_same_mass, tlog_t_same_mass};

        let args: [(R, R, R, R); 5] = [
            (0.03, 0.001, 2.12, 4.),
            (0.047, 0.001, 2.12, 4.),
            (0.047, 0.18, 2.12, 4.),
            (0.047, 0.18, 1.27, 4.),
            (0.047, 0.18, 1.27, 6.9),
        ];

        let res: [C; 5] = [
            1.142764958627341 + 0.0020362873111292638 * I,
            1.7889397625041612 + 0.003192663960369651 * I,
            1.8020841560466705 + 0.5746003840978272 * I,
            0.3920514512983847 + 0.34449821662631414 * I,
            0.39241163293773285 + 0.593714639016874 * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, en))| {
            assert_equal(tlog_t_same_mass(*q, *om, *p, *en), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, en))| {
            assert_equal(oneloop__tlog_t_same_mass(*q, *om, *p, *en), res[i])
        });
    }

    #[test]
    fn test_tlog_t_zero_mass() {
        use thermal::{ffi::oneloop__tlog_t_zero_mass, tlog_t_zero_mass};

        let args: [(R, R, R); 4] = [
            (0.03, 0.001, 2.12),
            (0.047, 0.001, 2.12),
            (0.047, 0.18, 2.12),
            (0.047, 0.18, 1.27),
        ];

        let res: [C; 4] = [
            1.142764956691111 + 1.5272154859334505e-05 * I,
            1.7889397550489738 + 3.751380169051372e-05 * I,
            1.8018695562823286 + 0.006752357819614953 * I,
            0.3914517606841969 + 0.004054041577866953 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p))| assert_equal(tlog_t_zero_mass(*q, *om, *p), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p))| {
            assert_equal(oneloop__tlog_t_zero_mass(*q, *om, *p), res[i])
        });
    }

    #[test]
    fn test_j_m_i() {
        use thermal::{ffi::oneloop__j_m_i, j_m_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            0.031125096655053638,
            0.025947316994062875,
            0.018328845146770915,
            0.002393477001341905,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(j_m_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__j_m_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_j_m_i_0() {
        use thermal::{ffi::oneloop__j_m_i, j_m_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.03397033162029322,
            0.026573487297514475,
            0.004111788232813987,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(j_m_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__j_m_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_j_0_i() {
        use thermal::{ffi::oneloop__j_0_i, j_0_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.03397033162029322,
            0.026573487297514475,
            0.004111788232813987,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(j_0_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__j_0_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_j_m_l_i() {
        use thermal::{ffi::oneloop__j_m_l_i, j_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            -0.004124075306794606,
            -0.013710562299662827,
            -0.01273488160797643,
            -0.0016629878205323555,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(j_m_l_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__j_m_l_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_j_m_l_i_0() {
        use thermal::{ffi::oneloop__j_m_l_i, j_m_l_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.0041613656234859185,
            -0.013775695815031503,
            -0.0021315510198907706,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(j_m_l_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__j_m_l_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_j_0_l_i() {
        use thermal::{ffi::oneloop__j_0_l_i, j_0_l_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.0041613656234859185,
            -0.013775695815031503,
            -0.0021315510198907706,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(j_0_l_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__j_0_l_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_j_m_t_i() {
        use thermal::{ffi::oneloop__j_m_t_i, j_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            0.0012709414467480234,
            0.004483696376574065,
            0.0031672244413620135,
            0.00041359282583188113,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(j_m_t_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__j_m_t_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_j_m_t_i_0() {
        use thermal::{ffi::oneloop__j_m_t_i, j_m_t_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.0013871218744953065,
            0.0045918986050105,
            0.0007105170066302567,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(j_m_t_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__j_m_t_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_j_0_t_i() {
        use thermal::{ffi::oneloop__j_0_t_i, j_0_t_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.0013871218744953065,
            0.0045918986050105,
            0.0007105170066302567,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(j_0_t_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__j_0_t_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_i_m_m_i() {
        use thermal::{ffi::oneloop__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            0.00021158661214306613 + 0. * I,
            0.00011022324250867698 + 0. * I,
            9.251825453605584e-05 + 0. * I,
            0.00015948317450630888 + 0. * I,
            2.990612245792503e-05 + 0. * I,
            1.7104766179985101e-07 + 0. * I,
            0.00029458934206651313 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(i_m_m_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__i_m_m_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_i_symm() {
        use thermal::{ffi::oneloop__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            0.00021158661214306613 + 0. * I,
            0.00011022324250867698 + 0. * I,
            9.251825453605584e-05 + 0. * I,
            0.00015948317450630888 + 0. * I,
            2.990612245792503e-05 + 0. * I,
            1.7104766179985101e-07 + 0. * I,
            0.00029458934206651313 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(i_m_m_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__i_m_m_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_i_m_0() {
        use thermal::{ffi::oneloop__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            0.00094234355160847 + 0. * I,
            0.0014541824000338444 + 0. * I,
            0.0013215297608047787 + 0. * I,
            0.002573892966401413 + 0. * I,
            0.000535813001876761 + 0. * I,
            0.01598483494461845 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_m_i(*q, *om, *p, *m, 0., *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_m_i(*q, *om, *p, *m, 0., *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_i_0_0() {
        use thermal::{ffi::oneloop__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            0.002400203493788911 + 0. * I,
            0.0037758168658999455 + 0. * I,
            0.003342530714848782 + 0. * I,
            0.008802715375634521 + 0. * I,
            0.2643407358846726 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(i_m_m_i(*q, *om, *p, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_m_m_i(*q, *om, *p, 0., 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_i_same_mass() {
        use thermal::{ffi::oneloop__i_m_m_i_same_mass, i_m_m_i_same_mass};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            9.096636827398976e-05 + 0. * I,
            4.0179611151967706e-05 + 0. * I,
            3.1811674811538385e-05 + 0. * I,
            4.782962588304627e-05 + 0. * I,
            9.668044889897316e-10 + 0. * I,
            0.0001631045339123352 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_m_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_m_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_0_i() {
        use thermal::{ffi::oneloop__i_m_0_i, i_m_0_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            0.00094234355160847 + 0. * I,
            0.0014541824000338444 + 0. * I,
            0.0013215297608047787 + 0. * I,
            0.002573892966401413 + 0. * I,
            0.000535813001876761 + 0. * I,
            0.01598483494461845 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_0_i(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_0_i(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_0_i_0_0() {
        use thermal::{ffi::oneloop__i_m_0_i, i_m_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            0.002400203493788911 + 0. * I,
            0.0037758168658999455 + 0. * I,
            0.003342530714848782 + 0. * I,
            0.008802715375634521 + 0. * I,
            0.2643407358846726 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(i_m_0_i(*q, *om, *p, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_m_0_i(*q, *om, *p, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_0_0_i() {
        use thermal::{ffi::oneloop__i_0_0_i, i_0_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            0.002400203493788911 + 0. * I,
            0.0037758168658999455 + 0. * I,
            0.003342530714848782 + 0. * I,
            0.008802715375634521 + 0. * I,
            0.2643407358846726 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, beta))| assert_equal(i_0_0_i(*q, *om, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_0_0_i(*q, *om, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_l_i() {
        use thermal::{ffi::oneloop__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            -0.000258329721743353 + 0. * I,
            -0.00010242087619035688 + 0. * I,
            -7.441517325762718e-05 + 0. * I,
            -8.942689983616806e-05 + 0. * I,
            -1.560319123456118e-05 + 0. * I,
            -5.972953834163304e-07 + 0. * I,
            -0.00017846136339926563 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(i_m_m_l_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__i_m_m_l_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_l_i_symm() {
        use thermal::{ffi::oneloop__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            -0.000258329721743353 + 0. * I,
            -0.00010242087619035688 + 0. * I,
            -7.441517325762718e-05 + 0. * I,
            -8.942689983616806e-05 + 0. * I,
            -1.560319123456118e-05 + 0. * I,
            -5.972953834163304e-07 + 0. * I,
            -0.00017846136339926563 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(i_m_m_l_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__i_m_m_l_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_l_i_m_0() {
        use thermal::{ffi::oneloop__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.0004569618607055666 + 0. * I,
            -0.00021774429666038944 + 0. * I,
            -0.00016214387969101272 + 0. * I,
            -0.00019054122091064876 + 0. * I,
            -3.9444985416917685e-05 + 0. * I,
            -0.00013405486365141085 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, beta))| {
                assert_equal(i_m_m_l_i(*q, *om, *p, *m1, 0., *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, beta))| {
                assert_equal(oneloop__i_m_m_l_i(*q, *om, *p, *m1, 0., *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_l_i_0_0() {
        use thermal::{ffi::oneloop__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.0009056407796552895 + 0. * I,
            -0.000456067961601411 + 0. * I,
            -0.0003451604129827074 + 0. * I,
            -0.0005621556195311336 + 0. * I,
            -0.01688122628159523 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(i_m_m_l_i(*q, *om, *p, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_m_m_l_i(*q, *om, *p, 0., 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_l_i_same_mass() {
        use thermal::{ffi::oneloop__i_m_m_l_i_same_mass, i_m_m_l_i_same_mass};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.00016337517007129169 + 0. * I,
            -6.208958860122221e-05 + 0. * I,
            -4.3657522785591495e-05 + 0. * I,
            -4.905839888881915e-05 + 0. * I,
            -9.5614472866962e-09 + 0. * I,
            -0.0016130618144456189 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_m_l_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_m_l_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_0_l_i() {
        use thermal::{ffi::oneloop__i_m_0_l_i, i_m_0_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.0004569618607055666 + 0. * I,
            -0.00021774429666038944 + 0. * I,
            -0.00016214387969101272 + 0. * I,
            -0.00019054122091064876 + 0. * I,
            -3.9444985416917685e-05 + 0. * I,
            -0.00013405486365141085 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_0_l_i(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_0_l_i(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_0_l_i_0_0() {
        use thermal::{ffi::oneloop__i_m_0_l_i, i_m_0_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.0009056407796552895 + 0. * I,
            -0.000456067961601411 + 0. * I,
            -0.0003451604129827074 + 0. * I,
            -0.0005621556195311336 + 0. * I,
            -0.01688122628159523 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(i_m_0_l_i(*q, *om, *p, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_m_0_l_i(*q, *om, *p, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_0_0_l_i() {
        use thermal::{ffi::oneloop__i_0_0_l_i, i_0_0_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.0009056407796552895 + 0. * I,
            -0.000456067961601411 + 0. * I,
            -0.0003451604129827074 + 0. * I,
            -0.0005621556195311336 + 0. * I,
            -0.01688122628159523 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, beta))| assert_equal(i_0_0_l_i(*q, *om, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_0_0_l_i(*q, *om, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_t_i() {
        use thermal::{ffi::oneloop__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            2.5928019321505628e-05 + 0. * I,
            4.445686894178127e-06 + 0. * I,
            3.7542561592958687e-06 + 0. * I,
            6.532283636844862e-06 + 0. * I,
            1.220736142629684e-06 + 0. * I,
            6.981446322597938e-09 + 0. * I,
            1.2018848342580335e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(i_m_m_t_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__i_m_m_t_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_t_i_symm() {
        use thermal::{ffi::oneloop__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            2.5928019321505628e-05 + 0. * I,
            4.445686894178127e-06 + 0. * I,
            3.7542561592958687e-06 + 0. * I,
            6.532283636844862e-06 + 0. * I,
            1.220736142629684e-06 + 0. * I,
            6.981446322597938e-09 + 0. * I,
            1.2018848342580335e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(i_m_m_t_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__i_m_m_t_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_t_i_m_0() {
        use thermal::{ffi::oneloop__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            0.0001168465635451649 + 0. * I,
            5.886366278185892e-05 + 0. * I,
            5.359851798598484e-05 + 0. * I,
            0.00010438792567407929 + 0. * I,
            2.1871697137444324e-05 + 0. * I,
            0.0006524944633346464 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, beta))| {
                assert_equal(i_m_m_t_i(*q, *om, *p, *m1, 0., *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, beta))| {
                assert_equal(oneloop__i_m_m_t_i(*q, *om, *p, *m1, 0., *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_m_t_i_0_0() {
        use thermal::{ffi::oneloop__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            0.00029181183135433764 + 0. * I,
            0.00015195021195586349 + 0. * I,
            0.0001349814205250193 + 0. * I,
            0.00035372481720296397 + 0. * I,
            0.01062216310422996 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(i_m_m_t_i(*q, *om, *p, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_m_m_t_i(*q, *om, *p, 0., 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_t_i_same_mass() {
        use thermal::{ffi::oneloop__i_m_m_t_i_same_mass, i_m_m_t_i_same_mass};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            1.111237552368281e-05 + 0. * I,
            1.6184614282070556e-06 + 0. * I,
            1.2922548805281137e-06 + 0. * I,
            1.9792396972389776e-06 + 0. * I,
            3.972600978190336e-11 + 0. * I,
            6.7019675471769906e-06 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_m_t_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_m_t_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_0_t_i() {
        use thermal::{ffi::oneloop__i_m_0_t_i, i_m_0_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            0.0001168465635451649 + 0. * I,
            5.886366278185892e-05 + 0. * I,
            5.359851798598484e-05 + 0. * I,
            0.00010438792567407929 + 0. * I,
            2.1871697137444324e-05 + 0. * I,
            0.0006524944633346464 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(i_m_0_t_i(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_0_t_i(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_i_m_0_t_i_0_0() {
        use thermal::{ffi::oneloop__i_m_0_t_i, i_m_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            0.00029181183135433764 + 0. * I,
            0.00015195021195586349 + 0. * I,
            0.0001349814205250193 + 0. * I,
            0.00035372481720296397 + 0. * I,
            0.01062216310422996 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(i_m_0_t_i(*q, *om, *p, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_m_0_t_i(*q, *om, *p, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_0_0_t_i() {
        use thermal::{ffi::oneloop__i_0_0_t_i, i_0_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            0.00029181183135433764 + 0. * I,
            0.00015195021195586349 + 0. * I,
            0.0001349814205250193 + 0. * I,
            0.00035372481720296397 + 0. * I,
            0.01062216310422996 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, beta))| assert_equal(i_0_0_t_i(*q, *om, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__i_0_0_t_i(*q, *om, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_d_j_m_i() {
        use thermal::{d_j_m_i, ffi::oneloop__d_j_m_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            -0.12704121083695363,
            -0.025026347409397066,
            -0.017678284285079972,
            -0.0023085233423436585,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(d_j_m_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__d_j_m_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_d_j_m_i_0() {
        use thermal::{d_j_m_i, ffi::oneloop__d_j_m_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.13865441477670704,
            -0.025630292532324916,
            -0.003965845131957935,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d_j_m_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d_j_m_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_d_j_0_i() {
        use thermal::{d_j_0_i, ffi::oneloop__d_j_0_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.13865441477670704,
            -0.025630292532324916,
            -0.003965845131957935,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d_j_0_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d_j_0_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_d2_j_m_i() {
        use thermal::{d2_j_m_i, ffi::oneloop__d2_j_m_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            1.071102648194892,
            0.05917388305807195,
            0.032849676717226224,
            0.006399344156501458,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(d2_j_m_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__d2_j_m_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_d2_j_m_i_0() {
        use thermal::{d2_j_m_i, ffi::oneloop__d2_j_m_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [1.259014323447687, 0.06163945774835798, 0.013570242883921217];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d2_j_m_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d2_j_m_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_d2_j_0_i() {
        use thermal::{d2_j_0_i, ffi::oneloop__d2_j_0_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [1.259014323447687, 0.06163945774835798, 0.013570242883921217];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d2_j_0_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d2_j_0_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_d_j_m_l_i() {
        use thermal::{d_j_m_l_i, ffi::oneloop__d_j_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            0.01683296043589635,
            0.013223921971125414,
            0.012282871921273563,
            0.001603962018260374,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(d_j_m_l_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__d_j_m_l_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_d_j_m_l_i_0() {
        use thermal::{d_j_m_l_i, ffi::oneloop__d_j_m_l_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.016985165810146613,
            0.013286743648757236,
            0.0020558941164069934,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d_j_m_l_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d_j_m_l_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_d_j_0_l_i() {
        use thermal::{d_j_0_l_i, ffi::oneloop__d_j_0_l_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.016985165810146613,
            0.013286743648757236,
            0.0020558941164069934,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d_j_0_l_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d_j_0_l_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_d2_j_m_l_i() {
        use thermal::{d2_j_m_l_i, ffi::oneloop__d2_j_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            -0.014879890048869512,
            -0.00624113239848815,
            -0.00514567109804881,
            -0.0021377409775935547,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(d2_j_m_l_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__d2_j_m_l_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_d2_j_m_l_i_0() {
        use thermal::{d2_j_m_l_i, ffi::oneloop__d2_j_m_l_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.015574839845634555,
            -0.006323602364423862,
            -0.003068968779066824,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d2_j_m_l_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d2_j_m_l_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_d2_j_0_l_i() {
        use thermal::{d2_j_0_l_i, ffi::oneloop__d2_j_0_l_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.015574839845634555,
            -0.006323602364423862,
            -0.003068968779066824,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d2_j_0_l_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d2_j_0_l_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_d_j_m_t_i() {
        use thermal::{d_j_m_t_i, ffi::oneloop__d_j_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            -0.015562548327526819,
            -0.012973658497031438,
            -0.009164422573385457,
            -0.0011967385006709525,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(d_j_m_t_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__d_j_m_t_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_d_j_m_t_i_0() {
        use thermal::{d_j_m_t_i, ffi::oneloop__d_j_m_t_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.01698516581014661,
            -0.013286743648757237,
            -0.0020558941164069934,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d_j_m_t_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d_j_m_t_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_d_j_0_t_i() {
        use thermal::{d_j_0_t_i, ffi::oneloop__d_j_0_t_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            -0.01698516581014661,
            -0.013286743648757237,
            -0.0020558941164069934,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d_j_0_t_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d_j_0_t_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_d2_j_m_t_i() {
        use thermal::{d2_j_m_t_i, ffi::oneloop__d2_j_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.35, 0.1, 1.2),
            (0.72, 0.1, 1.2),
            (0.72, 0.42, 1.2),
            (0.72, 0.42, 3.18),
        ];

        let res: [R; 4] = [
            0.06352060541847682,
            0.012513173704698533,
            0.008839142142539986,
            0.0011542616711718292,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(d2_j_m_t_i(*q, *m, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, m, beta))| assert_equal(oneloop__d2_j_m_t_i(*q, *m, *beta), res[i]))
    }

    #[test]
    fn test_d2_j_m_t_i_0() {
        use thermal::{d2_j_m_t_i, ffi::oneloop__d2_j_m_t_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.06932720738835352,
            0.012815146266162458,
            0.0019829225659789675,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d2_j_m_t_i(*q, 0., *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d2_j_m_t_i(*q, 0., *beta), res[i]))
    }

    #[test]
    fn test_d2_j_0_t_i() {
        use thermal::{d2_j_0_t_i, ffi::oneloop__d2_j_0_t_i};

        let args: [(R, R); 3] = [(0.35, 1.2), (0.72, 1.2), (0.72, 3.18)];

        let res: [R; 3] = [
            0.06932720738835352,
            0.012815146266162458,
            0.0019829225659789675,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(d2_j_0_t_i(*q, *beta), res[i]));

        args.iter()
            .enumerate()
            .for_each(|(i, (q, beta))| assert_equal(oneloop__d2_j_0_t_i(*q, *beta), res[i]))
    }

    #[test]
    fn test_d_i_m_m_i() {
        use thermal::{d_i_m_m_i, ffi::oneloop__d_i_m_m_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            -0.0002349261916872836 + 0. * I,
            -0.00023883399573377815 + 0. * I,
            -0.00016724617972627246 + 0. * I,
            -0.00021349292979744606 + 0. * I,
            -3.833673932041735e-06 - 0. * I,
            -1.926754710516342e-08 - 0. * I,
            0.0010328619019722522 - 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__d_i_m_m_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_i_symm() {
        use thermal::{d_i_m_m_i, ffi::oneloop__d_i_m_m_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            -0.0002349261916872836 + 0. * I,
            -0.00023883399573377815 + 0. * I,
            -0.00016724617972627246 + 0. * I,
            -0.00021349292979744606 + 0. * I,
            -3.833673932041735e-06 - 0. * I,
            -1.926754710516342e-08 - 0. * I,
            0.0010328619019722522 - 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__d_i_m_m_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_m_i_m_0() {
        use thermal::{d_i_m_m_i, ffi::oneloop__d_i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.0005364406452007649 + 0. * I,
            -0.0007386476734438202 + 0. * I,
            -0.0005553308703786081 + 0. * I,
            -0.001610573266508901 + 0. * I,
            -6.684519157056471e-05 + 0. * I,
            -0.0011357093446546923 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(d_i_m_m_i(*q, *om, *p, *m, 0., *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__d_i_m_m_i(*q, *om, *p, *m, 0., *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_m_i_0_0() {
        use thermal::{d_i_m_m_i, ffi::oneloop__d_i_m_m_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.004035658096018261 + 0. * I,
            -0.01656805554220818 + 0. * I,
            -0.014416047892326289 + 0. * I,
            -0.03880176518893119 + 0. * I,
            -1.1651958203779824 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(d_i_m_m_i(*q, *om, *p, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__d_i_m_m_i(*q, *om, *p, 0., 0., *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_i_same_mass() {
        use thermal::{d_i_m_m_i_same_mass, ffi::oneloop__d_i_m_m_i_same_mass};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.00014911215114362568 + 0. * I,
            -0.00017546859557974676 + 0. * I,
            -0.00013315954921477367 + 2.623437350273184e-21 * I,
            -0.0001818028795571539 + 1.0493749401092735e-20 * I,
            -3.823048863703415e-09 + 5.2587259851056895e-25 * I,
            -0.0006449665988725574 + 8.871721848019459e-20 * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(d_i_m_m_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__d_i_m_m_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_0_i() {
        use thermal::{d_i_m_0_i, ffi::oneloop__d_i_m_0_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -0.0005364406452007649 + 0. * I,
            -0.0007386476734438202 + 0. * I,
            -0.0005553308703786081 + 0. * I,
            -0.001610573266508901 + 0. * I,
            -6.684519157056471e-05 + 0. * I,
            -0.0011357093446546923 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(d_i_m_0_i(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__d_i_m_0_i(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_0_i_0_0() {
        use thermal::{d_i_m_0_i, ffi::oneloop__d_i_m_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.004035658096018261 + 0. * I,
            -0.01656805554220818 + 0. * I,
            -0.014416047892326289 + 0. * I,
            -0.03880176518893119 + 0. * I,
            -1.1651958203779824 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(d_i_m_0_i(*q, *om, *p, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__d_i_m_0_i(*q, *om, *p, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_0_0_i() {
        use thermal::{d_i_0_0_i, ffi::oneloop__d_i_0_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.004035658096018261 + 0. * I,
            -0.01656805554220818 + 0. * I,
            -0.014416047892326289 + 0. * I,
            -0.03880176518893119 + 0. * I,
            -1.1651958203779824 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, beta))| assert_equal(d_i_0_0_i(*q, *om, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__d_i_0_0_i(*q, *om, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_t_i() {
        use thermal::{d_i_m_m_t_i, ffi::oneloop__d_i_m_m_t_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            -3.0453698720331727e-05 + 0. * I,
            -1.259718018389033e-05 + 0. * I,
            -9.325325168374524e-06 + 0. * I,
            -1.1473523427992827e-05 + 0. * I,
            -7.79739481623871e-08 + 0. * I,
            -1.2851615018482306e-10 + 0. * I,
            6.569439555819919e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_t_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__d_i_m_m_t_i(*q, *om, *p, *m1, *m2, *beta), res[i])
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_t_i_symm() {
        use thermal::{d_i_m_m_t_i, ffi::oneloop__d_i_m_m_t_i};

        let args: [(R, R, R, R, R, R); 7] = [
            (0.62, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.21, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 2.16, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.8, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 3.2),
            (0.35, 0.75, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 7] = [
            -3.0453698720331727e-05 + 0. * I,
            -1.259718018389033e-05 + 0. * I,
            -9.325325168374524e-06 + 0. * I,
            -1.1473523427992827e-05 + 0. * I,
            -7.79739481623871e-08 + 0. * I,
            -1.2851615018482306e-10 + 0. * I,
            6.569439555819919e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_t_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m1, m2, beta))| {
                assert_equal(oneloop__d_i_m_m_t_i(*q, *om, *p, *m2, *m1, *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_m_t_i_m_0() {
        use thermal::{d_i_m_m_t_i, ffi::oneloop__d_i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -5.207781161469446e-05 + 0. * I,
            -2.4004003441611033e-05 + 0. * I,
            -1.773716674228576e-05 + 0. * I,
            -3.594307598059408e-05 + 0. * I,
            -1.3631787142386054e-06 + 0. * I,
            1.1960160876470227e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(d_i_m_m_t_i(*q, *om, *p, *m, 0., *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__d_i_m_m_t_i(*q, *om, *p, *m, 0., *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_m_t_i_0_0() {
        use thermal::{d_i_m_m_t_i, ffi::oneloop__d_i_m_m_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.0006000508734472277 + 0. * I,
            -0.0009439542164749866 + 0. * I,
            -0.0008356326787121956 + 0. * I,
            -0.0022006788439086302 + 0. * I,
            -0.06608518397116815 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(d_i_m_m_t_i(*q, *om, *p, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__d_i_m_m_t_i(*q, *om, *p, 0., 0., *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_t_i_same_mass() {
        use thermal::{d_i_m_m_t_i_same_mass, ffi::oneloop__d_i_m_m_t_i_same_mass};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -2.2741592068497444e-05 + 0. * I,
            -1.0044902787991927e-05 + 0. * I,
            -7.952918702884596e-06 + 0. * I,
            -1.1957406470761567e-05 + 0. * I,
            -2.417011222474329e-10 + 0. * I,
            -4.0776133478083794e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(d_i_m_m_t_i_same_mass(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(
                    oneloop__d_i_m_m_t_i_same_mass(*q, *om, *p, *m, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_d_i_m_0_t_i() {
        use thermal::{d_i_m_0_t_i, ffi::oneloop__d_i_m_0_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.21, 2.16, 1.2, 3.2),
            (0.35, 0.75, 2.16, 1.2, 3.2),
            (0.35, 0.75, 1.15, 1.2, 3.2),
            (0.35, 0.75, 1.15, 3.76, 3.2),
            (0.35, 0.75, 1.15, 3.76, 0.19),
        ];

        let res: [C; 6] = [
            -5.207781161469446e-05 + 0. * I,
            -2.4004003441611033e-05 + 0. * I,
            -1.773716674228576e-05 + 0. * I,
            -3.594307598059408e-05 + 0. * I,
            -1.3631787142386054e-06 + 0. * I,
            1.1960160876470227e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(d_i_m_0_t_i(*q, *om, *p, *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__d_i_m_0_t_i(*q, *om, *p, *m, *beta), res[i])
            })
    }

    #[test]
    fn test_d_i_m_0_t_i_0_0() {
        use thermal::{d_i_m_0_t_i, ffi::oneloop__d_i_m_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.0006000508734472277 + 0. * I,
            -0.0009439542164749866 + 0. * I,
            -0.0008356326787121956 + 0. * I,
            -0.0022006788439086302 + 0. * I,
            -0.06608518397116815 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(d_i_m_0_t_i(*q, *om, *p, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__d_i_m_0_t_i(*q, *om, *p, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_0_0_t_i() {
        use thermal::{d_i_0_0_t_i, ffi::oneloop__d_i_0_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 2.16, 3.2),
            (0.35, 0.21, 2.16, 3.2),
            (0.35, 0.75, 2.16, 3.2),
            (0.35, 0.75, 1.15, 3.2),
            (0.35, 0.75, 1.15, 0.19),
        ];

        let res: [C; 5] = [
            -0.0006000508734472277 + 0. * I,
            -0.0009439542164749866 + 0. * I,
            -0.0008356326787121956 + 0. * I,
            -0.0022006788439086302 + 0. * I,
            -0.06608518397116815 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(d_i_0_0_t_i(*q, *om, *p, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, beta))| {
            assert_equal(oneloop__d_i_0_0_t_i(*q, *om, *p, *beta), res[i])
        })
    }
}
