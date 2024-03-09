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

    pub fn tlog<T: Num>(q: R, p: R, delta: T) -> C {
        inlines::tlog(q, p, delta)
    }

    pub fn tlog_same_mass(q: R, p: R) -> C {
        inlines::tlog_same_mass(q, p)
    }

    pub fn tlog_zero_mass(q: R, p: R) -> C {
        inlines::tlog_zero_mass(q, p)
    }

    pub fn tlog_t<T: Num>(q: R, p: R, delta: T) -> C {
        inlines::tlog_t(q, p, delta)
    }

    pub fn tlog_t_same_mass(q: R, p: R) -> C {
        inlines::tlog_t_same_mass(q, p)
    }

    pub fn tlog_t_zero_mass(q: R, p: R) -> C {
        inlines::tlog_t_zero_mass(q, p)
    }

    pub fn i_m_m_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        inlines::i_m_m_i(q, p, m1, m2, beta)
    }

    pub fn i_m_m_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::i_m_m_i_same_mass(q, p, m, beta)
    }

    pub fn i_m_0_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::i_m_0_i(q, p, m, beta)
    }

    pub fn i_0_0_i(q: R, p: R, beta: R) -> C {
        inlines::i_0_0_i(q, p, beta)
    }

    pub fn i_m_m_l_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        inlines::i_m_m_l_i(q, p, m1, m2, beta)
    }

    pub fn i_m_m_l_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::i_m_m_l_i_same_mass(q, p, m, beta)
    }

    pub fn i_m_0_l_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::i_m_0_l_i(q, p, m, beta)
    }

    pub fn i_0_0_l_i(q: R, p: R, beta: R) -> C {
        inlines::i_0_0_l_i(q, p, beta)
    }

    pub fn i_m_m_t_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        inlines::i_m_m_t_i(q, p, m1, m2, beta)
    }

    pub fn i_m_m_t_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::i_m_m_t_i_same_mass(q, p, m, beta)
    }

    pub fn i_m_0_t_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::i_m_0_t_i(q, p, m, beta)
    }

    pub fn i_0_0_t_i(q: R, p: R, beta: R) -> C {
        inlines::i_0_0_t_i(q, p, beta)
    }

    pub fn d_i_m_m_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        inlines::d_i_m_m_i(q, p, m1, m2, beta)
    }

    pub fn d_i_m_m_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        inlines::d_i_m_m_i_same_mass(q, p, m, beta)
    }

    pub fn d_i_m_0_i<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        inlines::d_i_m_0_i(q, p, m, beta)
    }

    pub fn d_i_0_0_i(q: R, p: R, beta: R) -> R {
        inlines::d_i_0_0_i(q, p, beta)
    }

    pub fn d_i_m_m_l_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        inlines::d_i_m_m_l_i(q, p, m1, m2, beta)
    }

    pub fn d_i_m_m_l_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        inlines::d_i_m_m_l_i_same_mass(q, p, m, beta)
    }

    pub fn d_i_m_0_l_i<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        inlines::d_i_m_0_l_i(q, p, m, beta)
    }

    pub fn d_i_0_0_l_i(q: R, p: R, beta: R) -> R {
        inlines::d_i_0_0_l_i(q, p, beta)
    }

    pub fn d_i_m_m_t_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        inlines::d_i_m_m_t_i(q, p, m1, m2, beta)
    }

    pub fn d_i_m_m_t_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::d_i_m_m_t_i_same_mass(q, p, m, beta)
    }

    pub fn d_i_m_t_0_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, p, m, beta)
    }

    pub fn d_i_m_0_t_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, p, m, beta)
    }

    pub fn d_i_0_0_t_i(q: R, p: R, beta: R) -> C {
        inlines::d_i_0_0_t_i(q, p, beta)
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
    pub extern "C" fn oneloop__zero_matsubara__tlog(q: R, p: R, delta: R) -> C {
        inlines::tlog(q, p, delta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__tlog_same_mass(q: R, p: R) -> C {
        inlines::tlog_same_mass(q, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__tlog_zero_mass(q: R, p: R) -> C {
        inlines::tlog_zero_mass(q, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__tlog_t(q: R, p: R, delta: R) -> C {
        inlines::tlog_t(q, p, delta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__tlog_t_same_mass(q: R, p: R) -> C {
        inlines::tlog_t_same_mass(q, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__tlog_t_zero_mass(q: R, p: R) -> C {
        inlines::tlog_t_zero_mass(q, p)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_m_i(q: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_i(q, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_m_i_same_mass(q: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_m_i_same_mass(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_0_i(q: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_0_i(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_0_0_i(q: R, p: R, beta: R) -> C {
        inlines::i_0_0_i(q, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_m_l_i(q: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_l_i(q, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_m_l_i_same_mass(q: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_m_l_i_same_mass(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_0_l_i(q: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_0_l_i(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_0_0_l_i(q: R, p: R, beta: R) -> C {
        inlines::i_0_0_l_i(q, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_m_t_i(q: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::i_m_m_t_i(q, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_m_t_i_same_mass(q: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_m_t_i_same_mass(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_m_0_t_i(q: R, p: R, m: R, beta: R) -> C {
        inlines::i_m_0_t_i(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__i_0_0_t_i(q: R, p: R, beta: R) -> C {
        inlines::i_0_0_t_i(q, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_m_i(q: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_i(q, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_m_i_same_mass(q: R, p: R, m: R, beta: R) -> R {
        inlines::d_i_m_m_i_same_mass(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_0_i(q: R, p: R, m: R, beta: R) -> R {
        inlines::d_i_m_0_i(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_0_0_i(q: R, p: R, beta: R) -> R {
        inlines::d_i_0_0_i(q, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_m_l_i(q: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_l_i(q, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_m_l_i_same_mass(
        q: R,
        p: R,
        m: R,
        beta: R,
    ) -> R {
        inlines::d_i_m_m_l_i_same_mass(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_0_l_i(q: R, p: R, m: R, beta: R) -> R {
        inlines::d_i_m_0_l_i(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_0_0_l_i(q: R, p: R, beta: R) -> R {
        inlines::d_i_0_0_l_i(q, p, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_m_t_i(q: R, p: R, m1: R, m2: R, beta: R) -> C {
        inlines::d_i_m_m_t_i(q, p, m1, m2, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_m_t_i_same_mass(
        q: R,
        p: R,
        m: R,
        beta: R,
    ) -> C {
        inlines::d_i_m_m_t_i_same_mass(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_m_0_t_i(q: R, p: R, m: R, beta: R) -> C {
        inlines::d_i_m_0_t_i(q, p, m, beta)
    }

    #[no_mangle]
    pub extern "C" fn oneloop__zero_matsubara__d_i_0_0_t_i(q: R, p: R, beta: R) -> C {
        inlines::d_i_0_0_t_i(q, p, beta)
    }
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    //use num::traits::Inv;

    use crate::common::{inlines::energy, thermal::inlines::bose_distribution_zero_chempot};
    use crate::{Num, C, R};
    use std::f64::consts::PI;

    const PI2: R = PI * PI;

    #[inline(always)]
    pub fn r_0<T: Num>(p: R, delta: T) -> T {
        delta + p * p
    }

    // Note: this writing, while slowing down calculations, also avoids a branch cut
    // which does not exist in the Matsubara limit. It comes from
    // ln(A+i eta)+ln(A-i eta) = ln (A**2 + eta**2) -> ln(A**2) (as eta -> 0)
    // Note that the omega -> 0 limit of tlog in this module's parent module is not
    // well-defined precisely because of the branch cut (which cancels only when two
    // logs are summed), so this limit only makes sense when plugged into the integrals
    #[inline(always)]
    pub fn tlog<T: Num>(q: R, p: R, delta: T) -> C {
        let r_0 = r_0(p, delta);
        let q2p = q * 2. * p;
        let rat = (r_0 + q2p) / (r_0 - q2p);
        Into::<C>::into(rat * rat).ln() / 2.
    }

    #[inline(always)]
    pub fn tlog_same_mass(q: R, p: R) -> C {
        let rat = (p + 2. * q) / (p - 2. * q);
        (Into::<C>::into(rat * rat)).ln() / 2.
    }

    #[inline(always)]
    pub fn tlog_zero_mass(q: R, p: R) -> C {
        tlog_same_mass(q, p)
    }

    #[inline(always)]
    pub fn tlog_t<T: Num>(q: R, p: R, delta: T) -> C {
        let r0 = r_0(p, delta);
        let q2p = q * 2. * p;
        (r0 + q2p) * (r0 - q2p) * tlog(q, p, delta)
    }

    #[inline(always)]
    pub fn tlog_t_same_mass(q: R, p: R) -> C {
        p * p * (p + 2. * q) * (p - 2. * q) * tlog_same_mass(q, p)
    }

    #[inline(always)]
    pub fn tlog_t_zero_mass(q: R, p: R) -> C {
        tlog_t_same_mass(q, p)
    }

    #[inline(always)]
    pub fn i_m_m_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1 * tlog(q, p, delta);
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2 * tlog(q, p, -delta);
        (t1 + t2) * q / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        let t = bose_distribution_zero_chempot(en, beta) / en * tlog_same_mass(q, p);
        t * q / (p * 4. * PI2)
    }

    #[inline(always)]
    pub fn i_m_0_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        let delta = -m * m;
        let t1 = bose_distribution_zero_chempot(en, beta) * q / en * tlog(q, p, delta);
        let t2 = bose_distribution_zero_chempot(q, beta) * tlog(q, p, -delta);
        (t1 + t2) / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn i_0_0_i(q: R, p: R, beta: R) -> C {
        bose_distribution_zero_chempot(q, beta) * tlog_zero_mass(q, p) / (p * 4. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_l_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        -q * (bose_distribution_zero_chempot(en1, beta) * en1 * tlog(q, p, delta)
            + Into::<C>::into(bose_distribution_zero_chempot(en2, beta) * en2) * tlog(q, p, -delta))
            / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_l_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        -bose_distribution_zero_chempot(en, beta) * en * q * tlog_same_mass(q, p) / (p * 4. * PI2)
    }

    #[inline(always)]
    pub fn i_m_0_l_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        let delta = -m * m;
        -q * (bose_distribution_zero_chempot(en, beta) * en * tlog(q, p, delta)
            + Into::<C>::into(bose_distribution_zero_chempot(q, beta) * q) * tlog(q, p, -delta))
            / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn i_0_0_l_i(q: R, p: R, beta: R) -> C {
        -bose_distribution_zero_chempot(q, beta) * q * q * tlog_same_mass(q, p) / (p * 4. * PI2)
    }

    #[inline(always)]
    pub fn i_m_m_t_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let t1 = bose_distribution_zero_chempot(en1, beta) / en1
            * (tlog_t(q, p, delta) - 4. * qp * (p2 + delta));
        let t2 = bose_distribution_zero_chempot(en2, beta) / en2
            * (tlog_t(q, p, -delta) - 4. * qp * (p2 - delta));
        -(p3 * 64. * PI2).inv() * q * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_m_m_t_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let t = bose_distribution_zero_chempot(en, beta) / en
            * (-p2 * 4. * qp + tlog_t_same_mass(q, p));
        -(p3 * 32. * PI2).inv() * q * t
    }

    #[inline(always)]
    pub fn i_m_0_t_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        let delta = -Into::<C>::into(m * m);
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let t1 = bose_distribution_zero_chempot(en, beta) * q / en
            * (tlog_t(q, p, delta) - 4. * qp * (p2 + delta));
        let t2 = bose_distribution_zero_chempot(q, beta)
            * (tlog_t(q, p, -delta) - 4. * qp * (p2 - delta));
        -(p3 * 64. * PI2).inv() * (t1 + t2)
    }

    #[inline(always)]
    pub fn i_0_0_t_i(q: R, p: R, beta: R) -> C {
        let qp = q * p;
        let p2 = p * p;
        let p3 = p2 * p;
        let t = bose_distribution_zero_chempot(q, beta) * (-p2 * 4. * qp + tlog_t_zero_mass(q, p));
        -(p3 * 32. * PI2).inv() * t
    }

    #[inline(always)]
    pub fn d_i_m_m_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        let q2p = q * p * 2.;
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let en2 = energy(q, m2);
        let r01 = r_0(p, delta);
        let r02 = r_0(p, -delta);
        let r1pinv = (r01 + q2p).inv();
        let r1minv = (r01 - q2p).inv();
        let r2pinv = (r02 + q2p).inv();
        let r2minv = (r02 - q2p).inv();
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let t10 = b1 / en1 * (r1pinv + r1minv);
        let t1 = b1 / en1 * (r1pinv - r1minv);
        let t2 = b2 / en2 * (r2pinv - r2minv);
        (-t10 + (t2 - t1) * q / p) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        let rpinv = (p + 2. * q).inv();
        let rminv = (p - 2. * q).inv();
        -bose_distribution_zero_chempot(en, beta) * (en * p * 8. * PI2).inv() * (rpinv + rminv)
    }

    #[inline(always)]
    pub fn d_i_m_0_i<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        let q2p = q * p * 2.;
        let delta = -m * m;
        let en = energy(q, m);
        let r01 = r_0(p, delta);
        let r02 = r_0(p, -delta);
        let r1pinv = (r01 + q2p).inv();
        let r1minv = (r01 - q2p).inv();
        let r2pinv = (r02 + q2p).inv();
        let r2minv = (r02 - q2p).inv();
        let b1 = bose_distribution_zero_chempot(en, beta);
        let b2 = bose_distribution_zero_chempot(q, beta);
        let t10 = b1 / en * (r1pinv + r1minv);
        let t1 = b1 / en * (r1pinv - r1minv);
        let t2 = (r2pinv - r2minv) * b2 / q;
        (-t10 + (t2 - t1) * q / p) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_i(q: R, p: R, beta: R) -> R {
        let rpinv = (p + 2. * q).inv();
        let rminv = (p - 2. * q).inv();
        -bose_distribution_zero_chempot(q, beta) * (q * p * 8. * PI2).inv() * (rpinv + rminv)
    }

    #[inline(always)]
    pub fn d_i_m_m_l_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        let p2q = 2. * p * q;
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let b1en1 = b1 * en1;
        let r01 = r_0(p, delta);
        let r1pinv = (r01 + p2q).inv();
        let r1minv = (r01 - p2q).inv();
        let en2 = energy(q, m2);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let b2en2 = b2 * en2;
        let r02 = r_0(p, -delta);
        let r2pinv = (r02 + p2q).inv();
        let r2minv = (r02 - p2q).inv();
        let t1 = b1en1 * (r1pinv + r1minv);
        let t5 = -(b2en2 * (r2pinv - r2minv) - b1en1 * (r1pinv - r1minv)) * q / p;
        (t1 + t5) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_l_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        let en = energy(q, m);
        let r1pinv = (p + 2. * q).inv();
        let r1minv = (p - 2. * q).inv();
        bose_distribution_zero_chempot(en, beta) * en * (r1pinv + r1minv) / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_0_l_i<T: Num>(q: R, p: R, m: T, beta: R) -> T {
        let p2q = 2. * p * q;
        let delta = -m * m;
        let en = energy(q, m);
        let b1 = bose_distribution_zero_chempot(en, beta);
        let b1en1 = b1 * en;
        let r01 = r_0(p, delta);
        let r1pinv = (r01 + p2q).inv();
        let r1minv = (r01 - p2q).inv();
        let b2 = bose_distribution_zero_chempot(q, beta);
        let r02 = r_0(p, -delta);
        let r2pinv = (r02 + p2q).inv();
        let r2minv = (r02 - p2q).inv();
        let t1 = b1en1 * (r1pinv + r1minv);
        let t5 = -((r2pinv - r2minv) * b2 * q - b1en1 * (r1pinv - r1minv)) * q / p;
        (t1 + t5) / (8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_l_i(q: R, p: R, beta: R) -> R {
        let r1pinv = (p + 2. * q).inv();
        let r1minv = (p - 2. * q).inv();
        bose_distribution_zero_chempot(q, beta) * q * (r1pinv + r1minv) / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_t_i<T1: Num, T2: Num>(q: R, p: R, m1: T1, m2: T2, beta: R) -> C {
        let p2 = p * p;
        let delta = m2 * m2 - Into::<C>::into(m1 * m1);
        let en1 = energy(q, m1);
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let tlog1p = tlog(q, p, delta);
        let r01p = r_0(p, delta);
        let en2 = energy(q, m2);
        let b2 = bose_distribution_zero_chempot(en2, beta);
        let tlog2p = tlog(q, p, -delta);
        let r02p = r_0(p, -delta);
        let t1 = -b1 / en1 * tlog1p / 2.;
        let t2 = (b2 / en2 - Into::<C>::into(b1 / en1)) * q / p;
        let t3 = (b1 / en1 * (r01p * tlog1p) - b2 / en2 * (r02p * tlog2p)) / (p2 * 4.);
        (t1 + t2 + t3) * q / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_m_m_t_i_same_mass<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let en = energy(q, m);
        let b = bose_distribution_zero_chempot(en, beta);
        let tlogp = tlog_same_mass(q, p);
        -(en * p * 16. * PI2).inv() * b * tlogp * q
    }

    #[inline(always)]
    pub fn d_i_m_0_t_i<T: Num>(q: R, p: R, m: T, beta: R) -> C {
        let p2 = p * p;
        let delta = -m * m;
        let en1 = energy(q, m);
        let b1 = bose_distribution_zero_chempot(en1, beta);
        let tlog1p = tlog(q, p, delta);
        let r01p = r_0(p, delta);
        let b2 = bose_distribution_zero_chempot(q, beta);
        let tlog2p = tlog(q, p, -delta);
        let r02p = r_0(p, -delta);
        let t1 = -b1 / en1 * tlog1p / 2.;
        let t2 = (b2 / q - Into::<C>::into(b1 / en1)) * q / p;
        let t3 = (b1 / en1 * (r01p * tlog1p) - b2 / q * (r02p * tlog2p)) / (p2 * 4.);
        (t1 + t2 + t3) * q / (p * 8. * PI2)
    }

    #[inline(always)]
    pub fn d_i_0_0_t_i(q: R, p: R, beta: R) -> C {
        -bose_distribution_zero_chempot(q, beta) * tlog_same_mass(q, p) / (p * 16. * PI2)
    }
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::thermal::zero_matsubara;
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
        use zero_matsubara::{ffi::oneloop__zero_matsubara__tlog, tlog};

        let args: [(R, R, R); 4] = [
            (0.03, 2.12, 2.),
            (0.047, 2.12, 2.),
            (0.047, 1.27, 2.),
            (0.047, 1.27, 3.21),
        ];

        let res: [C; 4] = [
            0.039177220079546576 + 0. * I,
            0.06138906758007032 + 0. * I,
            0.06610948305432687 + 0. * I,
            0.04951559861177343 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, delta))| assert_equal(tlog(*q, *p, *delta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, delta))| {
            assert_equal(oneloop__zero_matsubara__tlog(*q, *p, *delta), res[i])
        });
    }

    #[test]
    fn test_tlog_same_mass() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__tlog_same_mass, tlog_same_mass};

        let args: [(R, R); 3] = [(0.03, 2.12), (0.047, 2.12), (0.047, 1.27)];

        let res: [C; 3] = [
            0.05661889399950811 + 0. * I,
            0.08873742845995349 + 0. * I,
            0.14830270994483524 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p))| assert_equal(tlog_same_mass(*q, *p), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p))| {
            assert_equal(oneloop__zero_matsubara__tlog_same_mass(*q, *p), res[i])
        });
    }

    #[test]
    fn test_tlog_zero_mass() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__tlog_zero_mass, tlog_zero_mass};

        let args: [(R, R); 3] = [(0.03, 2.12), (0.047, 2.12), (0.047, 1.27)];

        let res: [C; 3] = [
            0.05661889399950811 + 0. * I,
            0.08873742845995349 + 0. * I,
            0.14830270994483524 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p))| assert_equal(tlog_zero_mass(*q, *p), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p))| {
            assert_equal(oneloop__zero_matsubara__tlog_zero_mass(*q, *p), res[i])
        });
    }

    #[test]
    fn test_tlog_t() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__tlog_t, tlog_t};

        let args: [(R, R, R); 4] = [
            (0.03, 2.12, 2.),
            (0.047, 2.12, 2.),
            (0.047, 1.27, 2.),
            (0.047, 1.27, 3.21),
        ];

        let res: [C; 4] = [
            1.6517527941841421 + 0. * I,
            2.586782991823469 + 0. * I,
            0.861987985571848 + 0. * I,
            1.1510451919564078 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, delta))| assert_equal(tlog_t(*q, *p, *delta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, delta))| {
            assert_equal(oneloop__zero_matsubara__tlog_t(*q, *p, *delta), res[i])
        });
    }

    #[test]
    fn test_tlog_t_same_mass() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__tlog_t_same_mass, tlog_t_same_mass};

        let args: [(R, R); 3] = [(0.03, 2.12), (0.047, 2.12), (0.047, 1.27)];

        let res: [C; 3] = [
            1.1427647021550913 + 0. * I,
            1.788939355964949 + 0. * I,
            0.3836880037917354 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p))| assert_equal(tlog_t_same_mass(*q, *p), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p))| {
            assert_equal(oneloop__zero_matsubara__tlog_t_same_mass(*q, *p), res[i])
        });
    }

    #[test]
    fn test_tlog_t_zero_mass() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__tlog_t_zero_mass, tlog_t_zero_mass};

        let args: [(R, R); 3] = [(0.03, 2.12), (0.047, 2.12), (0.047, 1.27)];

        let res: [C; 3] = [
            1.1427647021550913 + 0. * I,
            1.788939355964949 + 0. * I,
            0.3836880037917354 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p))| assert_equal(tlog_t_zero_mass(*q, *p), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p))| {
            assert_equal(oneloop__zero_matsubara__tlog_t_zero_mass(*q, *p), res[i])
        });
    }

    #[test]
    fn test_i_m_m_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            0.00021895985426890845 + 0. * I,
            0.00011228270074367213 + 0. * I,
            0.0003186172300055374 - 0. * I,
            3.1269427573130054e-05 + 0. * I,
            1.9672643889409893e-07 + 0. * I,
            0.00029390621690845704 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(i_m_m_i(*q, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__i_m_m_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_i_symm() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            0.00021895985426890845 + 0. * I,
            0.00011228270074367213 + 0. * I,
            0.0003186172300055374 + 0. * I,
            3.1269427573130054e-05 + 0. * I,
            1.9672643889409893e-07 + 0. * I,
            0.00029390621690845704 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(i_m_m_i(*q, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__i_m_m_i(*q, *p, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_i_m_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_i, i_m_m_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.0009660137802164645 + 0. * I,
            0.0014674449940096626 + 0. * I,
            0.0031843131005532726 + 0. * I,
            0.0005559413779312267 + 0. * I,
            0.016572259609265443 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(i_m_m_i(*q, *p, *m, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_i(*q, *p, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_i_0_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_i, i_m_m_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.002444128404909455 + 0. * I,
            0.0038186975145049497 + 0. * I,
            0.01508023958599267 + 0. * I,
            0.4528513599919876 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_m_m_i(*q, *p, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_i(*q, *p, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_i_same_mass() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_i_same_mass, i_m_m_i_same_mass};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            9.462757686043952e-05 + 0. * I,
            4.119200029357628e-05 + 0. * I,
            0.00016266939999670026 + 0. * I,
            1.6303682659986545e-08 + 0. * I,
            0.0027505091169885637 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(i_m_m_i_same_mass(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_i_same_mass(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_0_i, i_m_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.0009660137802164645 + 0. * I,
            0.0014674449940096626 + 0. * I,
            0.0031843131005532726 + 0. * I,
            0.0005559413779312267 + 0. * I,
            0.016572259609265443 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(i_m_0_i(*q, *p, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(oneloop__zero_matsubara__i_m_0_i(*q, *p, *m, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_0_i_0_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_0_i, i_m_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.002444128404909455 + 0. * I,
            0.0038186975145049497 + 0. * I,
            0.01508023958599267 + 0. * I,
            0.4528513599919876 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_m_0_i(*q, *p, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__i_m_0_i(*q, *p, 0., *beta), res[i])
        })
    }

    #[test]
    fn test_i_0_0_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_0_0_i, i_0_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.002444128404909455 + 0. * I,
            0.0038186975145049497 + 0. * I,
            0.01508023958599267 + 0. * I,
            0.4528513599919876 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_0_0_i(*q, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__i_0_0_i(*q, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_l_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -0.0002737302427601515 + 0. * I,
            -0.00010586394672260547 + 0. * I,
            -0.0003141662984876339 + 0. * I,
            -2.3832631235225262e-05 + 0. * I,
            -9.82767094708467e-07 + 0. * I,
            0.0003507598523247177 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(i_m_m_l_i(*q, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__i_m_m_l_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_l_i_symm() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -0.0002737302427601515 + 0. * I,
            -0.00010586394672260547 + 0. * I,
            -0.0003141662984876339 + 0. * I,
            -2.3832631235225262e-05 + 0. * I,
            -9.82767094708467e-07 + 0. * I,
            0.0003507598523247177 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(i_m_m_l_i(*q, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__i_m_m_l_i(*q, *p, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_l_i_m_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.0004953119679504928 + 0. * I,
            -0.00022461796781799337 + 0. * I,
            -0.00036571886342540004 + 0. * I,
            -6.809256340417035e-05 + 0. * I,
            -0.00029996813256250835 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(i_m_m_l_i(*q, *p, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_l_i(*q, *p, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_l_i_0_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_l_i, i_m_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.0009395229588471947 + 0. * I,
            -0.0004677904455268564 + 0. * I,
            -0.0018473293492841018 + 0. * I,
            -0.055474291599018476 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_m_m_l_i(*q, *p, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_l_i(*q, *p, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_l_i_same_mass() {
        use zero_matsubara::{
            ffi::oneloop__zero_matsubara__i_m_m_l_i_same_mass, i_m_m_l_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.00017263855122418587 + 0. * I,
            -6.436250045871297e-05 + 0. * I,
            -0.0002541709374948442 + 0. * I,
            -2.3249214509967408e-07 + 0. * I,
            -0.03922253505916862 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(i_m_m_l_i_same_mass(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_l_i_same_mass(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_l_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_0_l_i, i_m_0_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.0004953119679504928 + 0. * I,
            -0.00022461796781799337 + 0. * I,
            -0.00036571886342540004 + 0. * I,
            -6.809256340417035e-05 + 0. * I,
            -0.00029996813256250835 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(i_m_0_l_i(*q, *p, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_0_l_i(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_l_i_0_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_0_l_i, i_m_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.0009395229588471947 + 0. * I,
            -0.0004677904455268564 + 0. * I,
            -0.0018473293492841018 + 0. * I,
            -0.055474291599018476 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_m_0_l_i(*q, *p, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_0_l_i(*q, *p, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_0_0_l_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_0_0_l_i, i_0_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.0009395229588471947 + 0. * I,
            -0.0004677904455268564 + 0. * I,
            -0.0018473293492841018 + 0. * I,
            -0.055474291599018476 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_0_0_l_i(*q, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__i_0_0_l_i(*q, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_i_m_m_t_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            2.661876938912228e-05 + 0. * I,
            4.524149407383567e-06 + 0. * I,
            1.691036349200717e-05 + 0. * I,
            1.2763318875993865e-06 + 0. * I,
            8.026749066028401e-09 + 0. * I,
            1.1996800492555959e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(i_m_m_t_i(*q, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__i_m_m_t_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_t_i_symm() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            2.661876938912228e-05 + 0. * I,
            4.524149407383567e-06 + 0. * I,
            1.691036349200717e-05 + 0. * I,
            1.2763318875993865e-06 + 0. * I,
            8.026749066028401e-09 + 0. * I,
            1.1996800492555959e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(i_m_m_t_i(*q, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__i_m_m_t_i(*q, *p, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_i_m_m_t_i_m_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00011882444004533821 + 0. * I,
            5.938278177055409e-05 + 0. * I,
            0.00012712441686931497 + 0. * I,
            2.269272158451893e-05 + 0. * I,
            0.0006764564522967755 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(i_m_m_t_i(*q, *p, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_t_i(*q, *p, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_t_i_0_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_m_t_i, i_m_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.0002963670345443 + 0. * I,
            0.00015361765012542159 + 0. * I,
            0.0005774993859039186 + 0. * I,
            0.01734199120709211 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_m_m_t_i(*q, *p, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_t_i(*q, *p, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_m_t_i_same_mass() {
        use zero_matsubara::{
            ffi::oneloop__zero_matsubara__i_m_m_t_i_same_mass, i_m_m_t_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            1.1474231175378945e-05 + 0. * I,
            1.657061934083353e-06 + 0. * I,
            6.229442050158883e-06 + 0. * I,
            6.243512691177725e-10 + 0. * I,
            0.00010533104046035392 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(i_m_m_t_i_same_mass(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_m_t_i_same_mass(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_t_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_0_t_i, i_m_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00011882444004533821 + 0. * I,
            5.938278177055409e-05 + 0. * I,
            0.00012712441686931497 + 0. * I,
            2.269272158451893e-05 + 0. * I,
            0.0006764564522967755 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(i_m_0_t_i(*q, *p, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_0_t_i(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_m_0_t_i_0_0() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_m_0_t_i, i_m_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.0002963670345443 + 0. * I,
            0.00015361765012542159 + 0. * I,
            0.0005774993859039186 + 0. * I,
            0.01734199120709211 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_m_0_t_i(*q, *p, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__i_m_0_t_i(*q, *p, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_i_0_0_t_i() {
        use zero_matsubara::{ffi::oneloop__zero_matsubara__i_0_0_t_i, i_0_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.0002963670345443 + 0. * I,
            0.00015361765012542159 + 0. * I,
            0.0005774993859039186 + 0. * I,
            0.01734199120709211 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(i_0_0_t_i(*q, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__i_0_0_t_i(*q, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_i() {
        use zero_matsubara::{d_i_m_m_i, ffi::oneloop__zero_matsubara__d_i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -0.00013512366005816255 + 0. * I,
            -0.000124935761700426 + 0. * I,
            0.000159728508866181 + 0. * I,
            -2.110957177353348e-06 + 0. * I,
            -1.3859086412622471e-08 + 0. * I,
            0.0007975402376583009 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_i(*q, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__d_i_m_m_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_i_symm() {
        use zero_matsubara::{d_i_m_m_i, ffi::oneloop__zero_matsubara__d_i_m_m_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -0.00013512366005816255 + 0. * I,
            -0.000124935761700426 + 0. * I,
            0.000159728508866181 + 0. * I,
            -2.110957177353348e-06 + 0. * I,
            -1.3859086412622471e-08 + 0. * I,
            0.0007975402376583009 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_i(*q, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__d_i_m_m_i(*q, *p, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_d_i_m_m_i_m_0() {
        use zero_matsubara::{d_i_m_m_i, ffi::oneloop__zero_matsubara__d_i_m_m_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -0.00035814806344659474 + 0. * I,
            -0.00038421028036718765 + 0. * I,
            -0.0014457519778219698 + 0. * I,
            -3.602192466474754e-05 + 0. * I,
            -0.0005714135841238835 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(d_i_m_m_i(*q, *p, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_i(*q, *p, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_i_0_0() {
        use zero_matsubara::{d_i_m_m_i, ffi::oneloop__zero_matsubara__d_i_m_m_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.002082531451785719 + 0. * I,
            -0.00839390284010033 + 0. * I,
            -0.04210159843738492 + 0. * I,
            -1.2642880109090309 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_m_m_i(*q, *p, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_i(*q, *p, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_i_same_mass() {
        use zero_matsubara::{
            d_i_m_m_i_same_mass, ffi::oneloop__zero_matsubara__d_i_m_m_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [R; 5] = [
            -8.06278854344546e-05,
            -9.054439293510987e-05,
            -0.0004541467473151286,
            -4.551725431729111e-08,
            -0.007678978154258283,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(d_i_m_m_i_same_mass(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_i_same_mass(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_i() {
        use zero_matsubara::{d_i_m_0_i, ffi::oneloop__zero_matsubara__d_i_m_0_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [R; 5] = [
            -0.00035814806344659474,
            -0.00038421028036718765,
            -0.0014457519778219698,
            -3.602192466474754e-05,
            -0.0005714135841238835,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(d_i_m_0_i(*q, *p, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_0_i(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_i_0_0() {
        use zero_matsubara::{d_i_m_0_i, ffi::oneloop__zero_matsubara__d_i_m_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [R; 4] = [
            -0.002082531451785719,
            -0.00839390284010033,
            -0.04210159843738492,
            -1.2642880109090309,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_m_0_i(*q, *p, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_0_i(*q, *p, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_0_0_i() {
        use zero_matsubara::{d_i_0_0_i, ffi::oneloop__zero_matsubara__d_i_0_0_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [R; 4] = [
            -0.002082531451785719,
            -0.00839390284010033,
            -0.04210159843738492,
            -1.2642880109090309,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_0_0_i(*q, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__d_i_0_0_i(*q, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_l_i() {
        use zero_matsubara::{d_i_m_m_l_i, ffi::oneloop__zero_matsubara__d_i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            0.00021884510543155269 + 0. * I,
            0.00018178617254283424 + 0. * I,
            -0.00034560348043263025 + 0. * I,
            1.5665621361903764e-06 + 0. * I,
            2.3898696881077647e-08 + 0. * I,
            -0.011805950974524424 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_l_i(*q, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__d_i_m_m_l_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_l_i_symm() {
        use zero_matsubara::{d_i_m_m_l_i, ffi::oneloop__zero_matsubara__d_i_m_m_l_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            0.00021884510543155269 + 0. * I,
            0.00018178617254283424 + 0. * I,
            -0.00034560348043263025 + 0. * I,
            1.5665621361903764e-06 + 0. * I,
            2.3898696881077647e-08 + 0. * I,
            -0.011805950974524424 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_l_i(*q, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m2, m1, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__d_i_m_m_l_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_d_i_m_m_l_i_m_0() {
        use zero_matsubara::{d_i_m_m_l_i, ffi::oneloop__zero_matsubara__d_i_m_m_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            0.00041385189120479195 + 0. * I,
            0.00024695011102484065 + 0. * I,
            0.0004882847601707977 + 0. * I,
            4.369914164071785e-06 + 0. * I,
            -0.007145776237643935 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(d_i_m_m_l_i(*q, *p, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_l_i(*q, *p, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_l_i_0_0() {
        use zero_matsubara::{d_i_m_m_l_i, ffi::oneloop__zero_matsubara__d_i_m_m_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            0.0008005250900664303 + 0. * I,
            0.0010282530979122905 + 0. * I,
            0.005157445808579652 + 0. * I,
            0.1548752813363563 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_m_m_l_i(*q, *p, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_l_i(*q, *p, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_l_i_same_mass() {
        use zero_matsubara::{
            d_i_m_m_l_i_same_mass, ffi::oneloop__zero_matsubara__d_i_m_m_l_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [R; 5] = [
            0.00014709751418661893,
            0.0001414756139611092,
            0.0007096042926798894,
            6.49080598290003e-07,
            0.10950299637753855,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(d_i_m_m_l_i_same_mass(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_l_i_same_mass(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_l_i() {
        use zero_matsubara::{d_i_m_0_l_i, ffi::oneloop__zero_matsubara__d_i_m_0_l_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [R; 5] = [
            0.00041385189120479195,
            0.00024695011102484065,
            0.0004882847601707977,
            4.369914164071785e-06,
            -0.007145776237643935,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(d_i_m_0_l_i(*q, *p, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_0_l_i(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_l_i_0_0() {
        use zero_matsubara::{d_i_m_0_l_i, ffi::oneloop__zero_matsubara__d_i_m_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [R; 4] = [
            0.0008005250900664303,
            0.0010282530979122905,
            0.005157445808579652,
            0.1548752813363563,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_m_0_l_i(*q, *p, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_0_l_i(*q, *p, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_0_0_l_i() {
        use zero_matsubara::{d_i_0_0_l_i, ffi::oneloop__zero_matsubara__d_i_0_0_l_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [R; 4] = [
            0.0008005250900664303,
            0.0010282530979122905,
            0.005157445808579652,
            0.1548752813363563,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_0_0_l_i(*q, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__d_i_0_0_l_i(*q, *p, *beta), res[i])
        })
    }

    #[test]
    fn test_d_i_m_m_t_i() {
        use zero_matsubara::{d_i_m_m_t_i, ffi::oneloop__zero_matsubara__d_i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -3.244028555061239e-05 + 0. * I,
            -1.303958804276184e-05 + 0. * I,
            -5.7830387408577956e-05 + 0. * I,
            -8.584123721940744e-08 + 0. * I,
            -1.7175969905878855e-10 + 0. * I,
            9.876174424143459e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_t_i(*q, *p, *m1, *m2, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__d_i_m_m_t_i(*q, *p, *m1, *m2, *beta),
                    res[i],
                )
            })
    }

    #[test]
    #[should_panic]
    fn test_d_i_m_m_t_i_symm() {
        use zero_matsubara::{d_i_m_m_t_i, ffi::oneloop__zero_matsubara__d_i_m_m_t_i};

        let args: [(R, R, R, R, R); 6] = [
            (0.62, 2.16, 1.2, 0.8, 3.2),
            (0.35, 2.16, 1.2, 0.8, 3.2),
            (0.35, 1.15, 1.2, 0.8, 3.2),
            (0.35, 1.15, 3.76, 0.8, 3.2),
            (0.35, 1.15, 3.76, 2.22, 3.2),
            (0.35, 1.15, 3.76, 2.22, 0.19),
        ];

        let res: [C; 6] = [
            -3.244028555061239e-05 + 0. * I,
            -1.303958804276184e-05 + 0. * I,
            -5.7830387408577956e-05 + 0. * I,
            -8.584123721940744e-08 + 0. * I,
            -1.7175969905878855e-10 + 0. * I,
            9.876174424143459e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(d_i_m_m_t_i(*q, *p, *m2, *m1, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m1, m2, beta))| {
                assert_equal(
                    oneloop__zero_matsubara__d_i_m_m_t_i(*q, *p, *m2, *m1, *beta),
                    res[i],
                )
            })
    }

    #[test]
    fn test_d_i_m_m_t_i_m_0() {
        use zero_matsubara::{d_i_m_m_t_i, ffi::oneloop__zero_matsubara__d_i_m_m_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -5.807951967299498e-05 + 0. * I,
            -2.4922506756341184e-05 + 0. * I,
            -5.711655480955709e-05 + 0. * I,
            -1.469056612934212e-06 + 0. * I,
            1.7453508136235085e-05 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(d_i_m_m_t_i(*q, *p, *m, 0., *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_t_i(*q, *p, *m, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_t_i_0_0() {
        use zero_matsubara::{d_i_m_m_t_i, ffi::oneloop__zero_matsubara__d_i_m_m_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.0006110321012273639 + 0. * I,
            -0.0009546743786262374 + 0. * I,
            -0.0037700598964981675 + 0. * I,
            -0.1132128399979969 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_m_m_t_i(*q, *p, 0., 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_t_i(*q, *p, 0., 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_m_t_i_same_mass() {
        use zero_matsubara::{
            d_i_m_m_t_i_same_mass, ffi::oneloop__zero_matsubara__d_i_m_m_t_i_same_mass,
        };

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -2.3656894215109883e-05 + 0. * I,
            -1.0298000073394072e-05 + 0. * I,
            -4.066734999917508e-05 + 0. * I,
            -4.075920664996636e-09 + 0. * I,
            -0.0006876272792471408 + 0. * I,
        ];

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(d_i_m_m_t_i_same_mass(*q, *p, *m, *beta), res[i])
        });

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_m_t_i_same_mass(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_t_i() {
        use zero_matsubara::{d_i_m_0_t_i, ffi::oneloop__zero_matsubara__d_i_m_0_t_i};

        let args: [(R, R, R, R); 5] = [
            (0.62, 2.16, 1.2, 3.2),
            (0.35, 2.16, 1.2, 3.2),
            (0.35, 1.15, 1.2, 3.2),
            (0.35, 1.15, 3.76, 3.2),
            (0.35, 1.15, 3.76, 0.19),
        ];

        let res: [C; 5] = [
            -5.807951967299498e-05 + 0. * I,
            -2.4922506756341184e-05 + 0. * I,
            -5.711655480955709e-05 + 0. * I,
            -1.469056612934212e-06 + 0. * I,
            1.7453508136235085e-05 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, m, beta))| assert_equal(d_i_m_0_t_i(*q, *p, *m, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, m, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_0_t_i(*q, *p, *m, *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_m_0_t_i_0_0() {
        use zero_matsubara::{d_i_m_0_t_i, ffi::oneloop__zero_matsubara__d_i_m_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.0006110321012273639 + 0. * I,
            -0.0009546743786262374 + 0. * I,
            -0.0037700598964981675 + 0. * I,
            -0.1132128399979969 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_m_0_t_i(*q, *p, 0., *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(
                oneloop__zero_matsubara__d_i_m_0_t_i(*q, *p, 0., *beta),
                res[i],
            )
        })
    }

    #[test]
    fn test_d_i_0_0_t_i() {
        use zero_matsubara::{d_i_0_0_t_i, ffi::oneloop__zero_matsubara__d_i_0_0_t_i};

        let args: [(R, R, R); 4] = [
            (0.62, 2.16, 3.2),
            (0.35, 2.16, 3.2),
            (0.35, 1.15, 3.2),
            (0.35, 1.15, 0.19),
        ];

        let res: [C; 4] = [
            -0.0006110321012273639 + 0. * I,
            -0.0009546743786262374 + 0. * I,
            -0.0037700598964981675 + 0. * I,
            -0.1132128399979969 + 0. * I,
        ];

        args.iter()
            .enumerate()
            .for_each(|(i, (q, p, beta))| assert_equal(d_i_0_0_t_i(*q, *p, *beta), res[i]));

        args.iter().enumerate().for_each(|(i, (q, p, beta))| {
            assert_equal(oneloop__zero_matsubara__d_i_0_0_t_i(*q, *p, *beta), res[i])
        })
    }
}
