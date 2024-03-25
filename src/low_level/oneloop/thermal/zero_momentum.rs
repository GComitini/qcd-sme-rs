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

    pub fn i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        inlines::i_m_m_i(q, om, m1, m2, beta)
    }

    pub fn i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::i_m_m_i_same_mass(q, om, m, beta)
    }

    pub fn i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::i_m_0_i(q, om, m, beta)
    }

    pub fn i_0_0_i<T: Num>(q: R, om: T, beta: R) -> C {
        inlines::i_0_0_i(q, om, beta)
    }

    pub fn i_m_m_l_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        inlines::i_m_m_l_i(q, om, m1, m2, beta)
    }

    pub fn i_m_m_l_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::i_m_m_l_i_same_mass(q, om, m, beta)
    }

    pub fn i_m_0_l_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::i_m_0_l_i(q, om, m, beta)
    }

    pub fn i_0_0_l_i<T: Num>(q: R, om: T, beta: R) -> C {
        inlines::i_0_0_l_i(q, om, beta)
    }

    pub fn i_m_m_t_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        inlines::i_m_m_t_i(q, om, m1, m2, beta)
    }

    pub fn i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::i_m_m_t_i_same_mass(q, om, m, beta)
    }

    pub fn i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::i_m_0_t_i(q, om, m, beta)
    }

    pub fn i_0_0_t_i<T: Num>(q: R, om: T, beta: R) -> C {
        inlines::i_0_0_t_i(q, om, beta)
    }

    pub fn d_i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        inlines::d_i_m_m_i(q, om, m1, m2, beta)
    }

    pub fn d_i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::d_i_m_m_i_same_mass(q, om, m, beta)
    }

    pub fn d_i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::d_i_m_0_i(q, om, m, beta)
    }

    pub fn d_i_0_0_i<T: Num>(q: R, om: T, beta: R) -> C {
        inlines::d_i_0_0_i(q, om, beta)
    }

    pub fn d_i_m_m_l_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        inlines::d_i_m_m_l_i(q, om, m1, m2, beta)
    }

    pub fn d_i_m_m_l_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::d_i_m_m_l_i_same_mass(q, om, m, beta)
    }

    pub fn d_i_m_0_l_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::d_i_m_0_l_i(q, om, m, beta)
    }

    pub fn d_i_0_0_l_i<T: Num>(q: R, om: T, beta: R) -> C {
        inlines::d_i_0_0_l_i(q, om, beta)
    }

    pub fn d_i_m_m_t_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        inlines::d_i_m_m_t_i(q, om, m1, m2, beta)
    }

    pub fn d_i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::d_i_m_m_t_i_same_mass(q, om, m, beta)
    }

    pub fn d_i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, om, m, beta)
    }

    pub fn d_i_0_0_t_i<T: Num>(q: R, om: T, beta: R) -> C {
        inlines::d_i_0_0_t_i(q, om, beta)
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
    pub extern "C" fn oneloop__zero_momentum__i_m_m_i(q: R, om: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_i(q, om, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_m_i_same_mass(q: R, om: R, m: R, beta: R) -> C {
        inlines::i_m_m_i_same_mass(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_0_i(q: R, om: R, m: R, beta: R) -> C {
        inlines::i_m_0_i(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_0_0_i(q: R, om: R, beta: R) -> C {
        inlines::i_0_0_i(q, om, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_m_l_i(q: R, om: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_l_i(q, om, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_m_l_i_same_mass(q: R, om: R, m: R, beta: R) -> C {
        inlines::i_m_m_l_i_same_mass(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_0_l_i(q: R, om: R, m: R, beta: R) -> C {
        inlines::i_m_0_l_i(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_0_0_l_i(q: R, om: R, beta: R) -> C {
        inlines::i_0_0_l_i(q, om, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_m_t_i(q: R, om: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_t_i(q, om, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_m_t_i_same_mass(q: R, om: R, m: R, beta: R) -> C {
        inlines::i_m_m_t_i_same_mass(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_m_0_t_i(q: R, om: R, m: R, beta: R) -> C {
        inlines::i_m_0_t_i(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__i_0_0_t_i(q: R, om: R, beta: R) -> C {
        inlines::i_0_0_t_i(q, om, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_m_i(q: R, om: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_i(q, om, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_m_i_same_mass(q: R, om: R, m: R, beta: R) -> C {
        inlines::d_i_m_m_i_same_mass(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_0_i(q: R, om: R, m: R, beta: R) -> C {
        inlines::d_i_m_0_i(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_0_0_i(q: R, om: R, beta: R) -> C {
        inlines::d_i_0_0_i(q, om, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_m_l_i(q: R, om: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_l_i(q, om, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_m_l_i_same_mass(
        q: R,
        om: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::d_i_m_m_l_i_same_mass(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_0_l_i(q: R, om: R, m: R, beta: R) -> C {
        inlines::d_i_m_0_l_i(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_0_0_l_i(q: R, om: R, beta: R) -> C {
        inlines::d_i_0_0_l_i(q, om, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_m_t_i(q: R, om: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_t_i(q, om, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_m_t_i_same_mass(
        q: R,
        om: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::d_i_m_m_t_i_same_mass(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_m_0_t_i(q: R, om: R, m: R, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, om, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_momentum__d_i_0_0_t_i(q: R, om: R, beta: R) -> C {
        inlines::d_i_0_0_t_i(q, om, beta)
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
    pub fn r_0<T1: Num, T2: Num, T3: Copy + Into<C>>(om: T1, en: T2, delta: T3) -> C {
        om * om + om * (en * 2. * I) + Into::<C>::into(delta)
    }

    #[inline(always)]
    pub fn r_0_same_mass<T1: Num, T2: Num>(om: T1, en: T2) -> C {
        om * om + om * (en * 2. * I)
    }

    #[inline(always)]
    pub fn i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (om2 + en1 * omi2 + delta).inv();
        let r01m_inv = (om2 - en1 * omi2 + delta).inv();
        let r02p_inv = (om2 + en2 * omi2 - delta).inv();
        let r02m_inv = (om2 - en2 * omi2 - delta).inv();
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1 * (r01p_inv + r01m_inv);
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2 * (r02p_inv + r02m_inv);
        (t1 + t2) * q * q / (4. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let en1 = energy(q, m);
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (om2 + en1 * omi2).inv();
        let r01m_inv = (om2 - en1 * omi2).inv();
        bose_distribution_zero_chempot(en1, beta) / en1 * (r01p_inv + r01m_inv) * q * q / (2. * PI2)
    }

    #[inline(always)]
    pub fn i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let en1 = energy(q, m);
        let delta = -m * m;
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (delta + (om2 + en1 * omi2)).inv();
        let r01m_inv = (delta + (om2 - en1 * omi2)).inv();
        let r02p_inv = (-delta + (om2 + q * omi2)).inv();
        let r02m_inv = (-delta + (om2 - q * omi2)).inv();
        let t1 = bose_distribution_zero_chempot(en1, beta) * q / en1 * (r01p_inv + r01m_inv);
        let t2 = bose_distribution_zero_chempot(q, beta) * (r02p_inv + r02m_inv);
        (t1 + t2) * q / (4. * PI2)
    }

    #[inline(always)]
    pub fn i_0_0_i<T: Num>(q: R, om: T, beta: R) -> C {
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (om2 + q * omi2).inv();
        let r01m_inv = (om2 - q * omi2).inv();
        bose_distribution_zero_chempot(q, beta) * (r01p_inv + r01m_inv) * q / (2. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_l_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let q4 = q * q * q * q;
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (om2 + en1 * omi2 + delta).inv();
        let r01m_inv = (om2 - en1 * omi2 + delta).inv();
        let r02p_inv = (om2 + en2 * omi2 - delta).inv();
        let r02m_inv = (om2 - en2 * omi2 - delta).inv();
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1 * (r01p_inv + r01m_inv);
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2 * (r02p_inv + r02m_inv);
        (12. * PI2).inv() * q4 * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_m_m_l_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let en1 = energy(q, m);
        let q4 = q * q * q * q;
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (om2 + en1 * omi2).inv();
        let r01m_inv = (om2 - en1 * omi2).inv();
        (en1 * 6. * PI2).inv()
            * bose_distribution_zero_chempot(en1, beta)
            * q4
            * (r01p_inv + r01m_inv)
    }

    #[inline(always)]
    pub fn i_m_0_l_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let en1 = energy(q, m);
        let delta = -m * m;
        let q3 = q * q * q;
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (delta + (om2 + en1 * omi2)).inv();
        let r01m_inv = (delta + (om2 - en1 * omi2)).inv();
        let r02p_inv = (-delta + (om2 + q * omi2)).inv();
        let r02m_inv = (-delta + (om2 - q * omi2)).inv();
        let t1 = bose_distribution_zero_chempot(en1, beta) * q / en1 * (r01p_inv + r01m_inv);
        let t2 = bose_distribution_zero_chempot(q, beta) * (r02p_inv + r02m_inv);
        (12. * PI2).inv() * q3 * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_0_0_l_i<T: Num>(q: R, om: T, beta: R) -> C {
        let q3 = q * q * q;
        let om2 = om * om;
        let omi2 = om * I * 2.;
        let r01p_inv = (om2 + q * omi2).inv();
        let r01m_inv = (om2 - q * omi2).inv();
        bose_distribution_zero_chempot(q, beta) * (r01p_inv + r01m_inv) * q3 / (6. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_t_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        i_m_m_l_i(q, om, m1, m2, beta)
    }

    #[inline(always)]
    pub fn i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        i_m_m_l_i_same_mass(q, om, m, beta)
    }

    #[inline(always)]
    pub fn i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        i_m_0_l_i(q, om, m, beta)
    }

    #[inline(always)]
    pub fn i_0_0_t_i<T: Num>(q: R, om: T, beta: R) -> C {
        i_0_0_l_i(q, om, beta)
    }

    #[inline(always)]
    pub fn d_i_m_m_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let r01_inv = r_0(om, en1, delta).inv();
        let r01_opp_inv = r_0(-om, en1, delta).inv();
        let r02_inv = r_0(om, en2, -delta).inv();
        let r02_opp_inv = r_0(-om, en2, -delta).inv();
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let t10 = b1 / en1 * (r01_inv + r01_opp_inv);
        let t1 = b1 / en1 * (r01_inv * r01_inv + r01_opp_inv * r01_opp_inv);
        let t2 = b2 / en2 * (r02_inv * r02_inv + r02_opp_inv * r02_opp_inv);
        (-t10 + (t1 - t2) * q * q * 2.) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let en = energy(q, m);
        let r0_inv = r_0_same_mass(om, en).inv();
        let r0_opp_inv = r_0_same_mass(-om, en).inv();
        -bose_distribution_zero_chempot(en, beta) * (en * 8. * PI2).inv() * (r0_inv + r0_opp_inv)
    }

    #[inline(always)]
    pub fn d_i_m_0_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let delta = -m * m;
        let en1 = energy(q, m);
        let r01_inv = r_0(om, en1, delta).inv();
        let r01_opp_inv = r_0(-om, en1, delta).inv();
        let r02_inv = r_0(om, q, -delta).inv();
        let r02_opp_inv = r_0(-om, q, -delta).inv();
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let b2 = bose_distribution_zero_chempot(q, beta);
        let t10 = b1 / en1 * (r01_inv + r01_opp_inv);
        let t1 = b1 * q / en1 * (r01_inv * r01_inv + r01_opp_inv * r01_opp_inv);
        let t2 = b2 * (r02_inv * r02_inv + r02_opp_inv * r02_opp_inv);
        (-t10 + (t1 - t2) * q * 2.) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_i<T: Num>(q: R, om: T, beta: R) -> C {
        let r0_inv = r_0_same_mass(om, q).inv();
        let r0_opp_inv = r_0_same_mass(-om, q).inv();
        let b1 = bose_distribution_zero_chempot(q, beta);
        -b1 * (r0_inv + r0_opp_inv) / (q * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_l_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        let q2 = q * q;
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let r01_inv = r_0(om, en1, delta).inv();
        let r01_opp_inv = r_0(-om, en1, delta).inv();
        let r02_inv = r_0(om, en2, -delta).inv();
        let r02_opp_inv = r_0(-om, en2, -delta).inv();
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let t10 = b1 / en1 * (r01_inv + r01_opp_inv);
        let t1 = b1 / en1 * (r01_inv * r01_inv + r01_opp_inv * r01_opp_inv);
        let t2 = b2 / en2 * (r02_inv * r02_inv + r02_opp_inv * r02_opp_inv);
        (-t10 + (t1 - t2) * q2 * 2. / 3.) * q2 / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_l_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let q2 = q * q;
        let en = energy(q, m);
        let r0_inv = r_0_same_mass(om, en).inv();
        let r0_opp_inv = r_0_same_mass(-om, en).inv();
        -bose_distribution_zero_chempot(en, beta)
            * (en * 8. * PI2).inv()
            * (r0_inv + r0_opp_inv)
            * q2
    }

    #[inline(always)]
    pub fn d_i_m_0_l_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        let q2 = q * q;
        let delta = -m * m;
        let en = energy(q, m);
        let r01_inv = r_0(om, en, delta).inv();
        let r01_opp_inv = r_0(-om, en, delta).inv();
        let r02_inv = r_0(om, q, -delta).inv();
        let r02_opp_inv = r_0(-om, q, -delta).inv();
        let b1 = bose_distribution_zero_chempot(en, beta);
        let b2 = bose_distribution_zero_chempot(q, beta);
        let t10 = b1 / en * (r01_inv + r01_opp_inv);
        let t1 = b1 * q / en * (r01_inv * r01_inv + r01_opp_inv * r01_opp_inv);
        let t2 = b2 * (r02_inv * r02_inv + r02_opp_inv * r02_opp_inv);
        (-t10 + (t1 - t2) * q * 2. / 3.) * q2 / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_l_i<T: Num>(q: R, om: T, beta: R) -> C {
        let r0_inv = r_0_same_mass(om, q).inv();
        let r0_opp_inv = r_0_same_mass(-om, q).inv();
        -bose_distribution_zero_chempot(q, beta) * (8. * PI2).inv() * (r0_inv + r0_opp_inv) * q
    }

    #[inline(always)]
    pub fn d_i_m_m_t_i<T1: Num, T2: Num, T3: Num>(q: R, om: T1, m1: T2, m2: T3, beta: R) -> C {
        d_i_m_m_l_i(q, om, m1, m2, beta)
    }

    #[inline(always)]
    pub fn d_i_m_m_t_i_same_mass<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        d_i_m_m_l_i_same_mass(q, om, m, beta)
    }

    #[inline(always)]
    pub fn d_i_m_0_t_i<T1: Num, T2: Num>(q: R, om: T1, m: T2, beta: R) -> C {
        d_i_m_0_l_i(q, om, m, beta)
    }

    #[inline(always)]
    pub fn d_i_0_0_t_i<T: Num>(q: R, om: T, beta: R) -> C {
        d_i_0_0_l_i(q, om, beta)
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::thermal::zero_momentum;
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
    fn test_i_m_m_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            0.0005779677841385296 + 0. * I,
            0.0003787155935185114 + 0. * I,
            0.00017043494005206962 + 0. * I,
            3.2644108168296024e-05 + 0. * I,
            1.8949064200227601e-07 + 0. * I,
            0.0003473435183059436 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(i_m_m_i(*q, *om, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__i_m_m_i(*q, *om, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_i_symm() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            0.0005779677841385296 + 0. * I,
            0.0003787155935185114 + 0. * I,
            0.00017043494005206962 + 0. * I,
            3.2644108168296024e-05 + 0. * I,
            1.8949064200227601e-07 + 0. * I,
            0.0003473435183059436 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(i_m_m_i(*q, *om, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__i_m_m_i(*q, *om, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_i_m_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.003154494690178643 + 0. * I,
            0.005671747675363454 + 0. * I,
            0.0039934561672281006 + 0. * I,
            0.0005834109847153382 + 0. * I,
            0.01742126492895578 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(i_m_m_i(*q, *om, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_i(*q, *om, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_i_0_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_i, i_m_m_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.006332534568446089 + 0. * I,
            0.032155578623048575 + 0. * I,
            0.016317619517881465 + 0. * I,
            0.4900092036580204 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_m_m_i(*q, *om, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_i(*q, *om, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_i_same_mass() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_i_same_mass, i_m_m_i_same_mass};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            5.2820043017623486e-05 + 0. * I,
            2.9433566657429e-05 + 0. * I,
            2.7193807251159485e-05 + 0. * I,
            3.2233838864676903e-10 + 0. * I,
            5.4380025373302544e-05 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(i_m_m_i_same_mass(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_i_same_mass(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_0_i, i_m_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.003154494690178643 + 0. * I,
            0.005671747675363454 + 0. * I,
            0.0039934561672281006 + 0. * I,
            0.0005834109847153382 + 0. * I,
            0.01742126492895578 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m, beta))| assert_equal(i_m_0_i(*q, *om, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(oneloop__zero_momentum__i_m_0_i(*q, *om, *m, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_0_i_0_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_0_i, i_m_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.006332534568446089 + 0. * I,
            0.032155578623048575 + 0. * I,
            0.016317619517881465 + 0. * I,
            0.4900092036580204 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_m_0_i(*q, *om, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__i_m_0_i(*q, *om, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_0_0_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_0_0_i, i_0_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.006332534568446089 + 0. * I,
            0.032155578623048575 + 0. * I,
            0.016317619517881465 + 0. * I,
            0.4900092036580204 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_0_0_i(*q, *om, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__i_0_0_i(*q, *om, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_l_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            7.405693874095026e-05 + 0. * I,
            1.546422006867255e-05 + 0. * I,
            6.959426718792841e-06 + 0. * I,
            1.3329677502054206e-06 + 0. * I,
            7.737534548426272e-09 + 0. * I,
            1.4183193664159363e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(i_m_m_l_i(*q, *om, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__i_m_m_l_i(*q, *om, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_l_i_symm() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            7.405693874095026e-05 + 0. * I,
            1.546422006867255e-05 + 0. * I,
            6.959426718792841e-06 + 0. * I,
            1.3329677502054206e-06 + 0. * I,
            7.737534548426272e-09 + 0. * I,
            1.4183193664159363e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(i_m_m_l_i(*q, *om, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__i_m_m_l_i(*q, *om, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_l_i_m_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00040419591963489015 + 0. * I,
            0.00023159636341067437 + 0. * I,
            0.00016306612682848078 + 0. * I,
            2.3822615209209636e-05 + 0. * I,
            0.0007113683179323609 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(i_m_m_l_i(*q, *om, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_l_i(*q, *om, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_l_i_0_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.0008114087627035585 + 0. * I,
            0.00131301946044115 + 0. * I,
            0.0006663027969801598 + 0. * I,
            0.020008709149369166 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_m_m_l_i(*q, *om, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_l_i(*q, *om, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_l_i_same_mass() {
        use zero_momentum::{
            ffi::oneloop__zero_momentum__i_m_m_l_i_same_mass, i_m_m_l_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            6.768008178658155e-06 + 0. * I,
            1.2018706385116842e-06 + 0. * I,
            1.110413796089012e-06 + 0. * I,
            1.3162150869743068e-11 + 0. * I,
            2.2205177027431874e-06 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(i_m_m_l_i_same_mass(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_l_i_same_mass(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_l_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_0_l_i, i_m_0_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00040419591963489015 + 0. * I,
            0.00023159636341067437 + 0. * I,
            0.00016306612682848078 + 0. * I,
            2.3822615209209636e-05 + 0. * I,
            0.0007113683179323609 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m, beta))| assert_equal(i_m_0_l_i(*q, *om, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_0_l_i(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_l_i_0_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_0_l_i, i_m_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.0008114087627035585 + 0. * I,
            0.00131301946044115 + 0. * I,
            0.0006663027969801598 + 0. * I,
            0.020008709149369166 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_m_0_l_i(*q, *om, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_0_l_i(*q, *om, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_0_0_l_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_0_0_l_i, i_0_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.0008114087627035585 + 0. * I,
            0.00131301946044115 + 0. * I,
            0.0006663027969801598 + 0. * I,
            0.020008709149369166 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_0_0_l_i(*q, *om, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__i_0_0_l_i(*q, *om, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_t_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            7.405693874095026e-05 + 0. * I,
            1.546422006867255e-05 + 0. * I,
            6.959426718792841e-06 + 0. * I,
            1.3329677502054206e-06 + 0. * I,
            7.737534548426272e-09 + 0. * I,
            1.4183193664159363e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(i_m_m_t_i(*q, *om, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__i_m_m_t_i(*q, *om, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_t_i_symm() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            7.405693874095026e-05 + 0. * I,
            1.546422006867255e-05 + 0. * I,
            6.959426718792841e-06 + 0. * I,
            1.3329677502054206e-06 + 0. * I,
            7.737534548426272e-09 + 0. * I,
            1.4183193664159363e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(i_m_m_t_i(*q, *om, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__i_m_m_t_i(*q, *om, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_t_i_m_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00040419591963489015 + 0. * I,
            0.00023159636341067437 + 0. * I,
            0.00016306612682848078 + 0. * I,
            2.3822615209209636e-05 + 0. * I,
            0.0007113683179323609 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(i_m_m_t_i(*q, *om, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_t_i(*q, *om, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_t_i_0_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.0008114087627035585 + 0. * I,
            0.00131301946044115 + 0. * I,
            0.0006663027969801598 + 0. * I,
            0.020008709149369166 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_m_m_t_i(*q, *om, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_t_i(*q, *om, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_t_i_same_mass() {
        use zero_momentum::{
            ffi::oneloop__zero_momentum__i_m_m_t_i_same_mass, i_m_m_t_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            6.768008178658155e-06 + 0. * I,
            1.2018706385116842e-06 + 0. * I,
            1.110413796089012e-06 + 0. * I,
            1.3162150869743068e-11 + 0. * I,
            2.2205177027431874e-06 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(i_m_m_t_i_same_mass(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_m_t_i_same_mass(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_t_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_0_t_i, i_m_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00040419591963489015 + 0. * I,
            0.00023159636341067437 + 0. * I,
            0.00016306612682848078 + 0. * I,
            2.3822615209209636e-05 + 0. * I,
            0.0007113683179323609 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m, beta))| assert_equal(i_m_0_t_i(*q, *om, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_0_t_i(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_t_i_0_0() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_m_0_t_i, i_m_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.0008114087627035585 + 0. * I,
            0.00131301946044115 + 0. * I,
            0.0006663027969801598 + 0. * I,
            0.020008709149369166 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_m_0_t_i(*q, *om, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__i_m_0_t_i(*q, *om, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_0_0_t_i() {
        use zero_momentum::{ffi::oneloop__zero_momentum__i_0_0_t_i, i_0_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            0.0008114087627035585 + 0. * I,
            0.00131301946044115 + 0. * I,
            0.0006663027969801598 + 0. * I,
            0.020008709149369166 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(i_0_0_t_i(*q, *om, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__i_0_0_t_i(*q, *om, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_i() {
        use zero_momentum::{d_i_m_m_i, ffi::oneloop__zero_momentum__d_i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -0.00024932117296966094 + 0. * I,
            2.598618215220097e-06 + 0. * I,
            -5.083053176132579e-06 + 0. * I,
            -2.2793554053238165e-06 + 0. * I,
            -1.2226188761266224e-08 + 0. * I,
            0.0004849808790824091 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(d_i_m_m_i(*q, *om, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__d_i_m_m_i(*q, *om, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_i_symm() {
        use zero_momentum::{d_i_m_m_i, ffi::oneloop__zero_momentum__d_i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -0.00024932117296966094 + 0. * I,
            2.598618215220097e-06 + 0. * I,
            -5.083053176132579e-06 + 0. * I,
            -2.2793554053238165e-06 + 0. * I,
            -1.2226188761266224e-08 + 0. * I,
            0.0004849808790824091 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(d_i_m_m_i(*q, *om, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__d_i_m_m_i(*q, *om, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_d_i_m_m_i_m_0() {
        use zero_momentum::{d_i_m_m_i, ffi::oneloop__zero_momentum__d_i_m_m_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.0018578359639273342 + 0. * I,
            -0.0035172171691799023 + 0. * I,
            -0.0016824090211125752 + 0. * I,
            -3.9584070600586926e-05 + 0. * I,
            -0.0007825833677795731 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_m_i(*q, *om, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_i(*q, *om, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_i_0_0() {
        use zero_momentum::{d_i_m_m_i, ffi::oneloop__zero_momentum__d_i_m_m_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.004118453803619982 + 0. * I,
            -0.06562362984295625 + 0. * I,
            -0.03330126432220708 + 0. * I,
            -1.000018782975552 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_m_m_i(*q, *om, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_i(*q, *om, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_i_same_mass() {
        use zero_momentum::{
            d_i_m_m_i_same_mass, ffi::oneloop__zero_momentum__d_i_m_m_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -3.4352265229984086e-05 + 0. * I,
            -6.006850338250823e-05 + 0. * I,
            -5.5497565818692824e-05 + 0. * I,
            -6.578334462178961e-10 + 0. * I,
            -0.0001109796436189848 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_m_i_same_mass(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_i_same_mass(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_i() {
        use zero_momentum::{d_i_m_0_i, ffi::oneloop__zero_momentum__d_i_m_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.0018578359639273342 + 0. * I,
            -0.0035172171691799023 + 0. * I,
            -0.0016824090211125752 + 0. * I,
            -3.9584070600586926e-05 + 0. * I,
            -0.0007825833677795731 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m, beta))| assert_equal(d_i_m_0_i(*q, *om, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_0_i(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_i_0_0() {
        use zero_momentum::{d_i_m_0_i, ffi::oneloop__zero_momentum__d_i_m_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.004118453803619982 + 0. * I,
            -0.06562362984295625 + 0. * I,
            -0.03330126432220708 + 0. * I,
            -1.000018782975552 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_m_0_i(*q, *om, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_0_i(*q, *om, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_0_0_i() {
        use zero_momentum::{d_i_0_0_i, ffi::oneloop__zero_momentum__d_i_0_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.004118453803619982 + 0. * I,
            -0.06562362984295625 + 0. * I,
            -0.03330126432220708 + 0. * I,
            -1.000018782975552 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_0_0_i(*q, *om, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__d_i_0_0_i(*q, *om, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_l_i() {
        use zero_momentum::{d_i_m_m_l_i, ffi::oneloop__zero_momentum__d_i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            2.2749626510847174e-05 + 0. * I,
            2.766115906785751e-05 + 0. * I,
            1.8453696456422319e-06 + 0. * I,
            -9.287293377729397e-08 + 0. * I,
            -2.488106457198816e-10 + 0. * I,
            6.205134590847087e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(d_i_m_m_l_i(*q, *om, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__d_i_m_m_l_i(*q, *om, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_l_i_symm() {
        use zero_momentum::{d_i_m_m_l_i, ffi::oneloop__zero_momentum__d_i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.21, 1.2, 0.8, 3.2),
            (0.35, 0.75, 1.2, 0.8, 3.2),
            (0.35, 0.75, 3.76, 0.8, 3.2),
            (0.35, 0.75, 3.76, 2.22, 3.2),
            (0.35, 0.75, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            2.2749626510847174e-05 + 0. * I,
            2.766115906785751e-05 + 0. * I,
            1.8453696456422319e-06 + 0. * I,
            -9.287293377729397e-08 + 0. * I,
            -2.488106457198816e-10 + 0. * I,
            6.205134590847087e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(d_i_m_m_l_i(*q, *om, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_momentum__d_i_m_m_l_i(*q, *om, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_d_i_m_m_l_i_m_0() {
        use zero_momentum::{d_i_m_m_l_i, ffi::oneloop__zero_momentum__d_i_m_m_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.00019831294139692566 + 0. * I,
            -0.0001242415111677944 + 0. * I,
            -6.237632304143158e-05 + 0. * I,
            -1.6161553920999904e-06 + 0. * I,
            7.997970344844774e-07 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_m_l_i(*q, *om, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_l_i(*q, *om, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_l_i_0_0() {
        use zero_momentum::{d_i_m_m_l_i, ffi::oneloop__zero_momentum__d_i_m_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.001583133642111522 + 0. * I,
            -0.008038894655762144 + 0. * I,
            -0.004079404879470366 + 0. * I,
            -0.1225023009145051 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(d_i_m_m_l_i(*q, *om, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_l_i(*q, *om, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_l_i_same_mass() {
        use zero_momentum::{
            d_i_m_m_l_i_same_mass, ffi::oneloop__zero_momentum__d_i_m_m_l_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -1.3205010754405883e-05 + 0. * I,
            -7.358391664357257e-06 + 0. * I,
            -6.7984518127898696e-06 + 0. * I,
            -8.058459716169227e-11 + 0. * I,
            -1.3595006343325636e-05 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_m_l_i_same_mass(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_l_i_same_mass(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_l_i() {
        use zero_momentum::{d_i_m_0_l_i, ffi::oneloop__zero_momentum__d_i_m_0_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.00019831294139692566 + 0. * I,
            -0.0001242415111677944 + 0. * I,
            -6.237632304143158e-05 + 0. * I,
            -1.6161553920999904e-06 + 0. * I,
            7.997970344844774e-07 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_0_l_i(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_0_l_i(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_l_i_0_0() {
        use zero_momentum::{d_i_m_0_l_i, ffi::oneloop__zero_momentum__d_i_m_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.001583133642111522 + 0. * I,
            -0.008038894655762144 + 0. * I,
            -0.004079404879470366 + 0. * I,
            -0.1225023009145051 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_m_0_l_i(*q, *om, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_0_l_i(*q, *om, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_0_0_l_i() {
        use zero_momentum::{d_i_0_0_l_i, ffi::oneloop__zero_momentum__d_i_0_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.001583133642111522 + 0. * I,
            -0.008038894655762144 + 0. * I,
            -0.004079404879470366 + 0. * I,
            -0.1225023009145051 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_0_0_l_i(*q, *om, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__d_i_0_0_l_i(*q, *om, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_t_i_0_0() {
        use zero_momentum::{d_i_m_m_t_i, ffi::oneloop__zero_momentum__d_i_m_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.001583133642111522 + 0. * I,
            -0.008038894655762144 + 0. * I,
            -0.004079404879470366 + 0. * I,
            -0.1225023009145051 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(d_i_m_m_t_i(*q, *om, 0., 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_t_i(*q, *om, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_t_i_same_mass() {
        use zero_momentum::{
            d_i_m_m_t_i_same_mass, ffi::oneloop__zero_momentum__d_i_m_m_t_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -1.3205010754405883e-05 + 0. * I,
            -7.358391664357257e-06 + 0. * I,
            -6.7984518127898696e-06 + 0. * I,
            -8.058459716169227e-11 + 0. * I,
            -1.3595006343325636e-05 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_m_t_i_same_mass(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_m_t_i_same_mass(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_t_i() {
        use zero_momentum::{d_i_m_0_t_i, ffi::oneloop__zero_momentum__d_i_m_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 0.21, 1.2, 3.2),
            (0.35, 0.21, 1.2, 3.2),
            (0.35, 0.75, 1.2, 3.2),
            (0.35, 0.75, 3.76, 3.2),
            (0.35, 0.75, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.00019831294139692566 + 0. * I,
            -0.0001242415111677944 + 0. * I,
            -6.237632304143158e-05 + 0. * I,
            -1.6161553920999904e-06 + 0. * I,
            7.997970344844774e-07 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(d_i_m_0_t_i(*q, *om, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, om, m, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_0_t_i(*q, *om, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_t_i_0_0() {
        use zero_momentum::{d_i_m_0_t_i, ffi::oneloop__zero_momentum__d_i_m_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.001583133642111522 + 0. * I,
            -0.008038894655762144 + 0. * I,
            -0.004079404879470366 + 0. * I,
            -0.1225023009145051 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_m_0_t_i(*q, *om, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(
                oneloop__zero_momentum__d_i_m_0_t_i(*q, *om, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_0_0_t_i() {
        use zero_momentum::{d_i_0_0_t_i, ffi::oneloop__zero_momentum__d_i_0_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 0.21, 3.2),
            (0.35, 0.21, 3.2),
            (0.35, 0.75, 3.2),
            (0.35, 0.75, 0.19),
        ];

        let res: [C; 4] = [
            -0.001583133642111522 + 0. * I,
            -0.008038894655762144 + 0. * I,
            -0.004079404879470366 + 0. * I,
            -0.1225023009145051 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, beta))| assert_equal(d_i_0_0_t_i(*q, *om, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, om, beta))| {
            assert_equal(oneloop__zero_momentum__d_i_0_0_t_i(*q, *om, *beta), res[i])
        })
    }
}
