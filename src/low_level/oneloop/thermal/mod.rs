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
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
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
        ((r_0(om, p, en, delta) + q * 2. * p) / (r_0(om, p, en, delta) - q * 2. * p)).ln()
    }

    #[inline(always)]
    pub fn tlog_same_mass<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2) -> C {
        ((r_0_same_mass(om, p, en) + q * 2. * p) / (r_0_same_mass(om, p, en) - q * 2. * p)).ln()
    }

    #[inline(always)]
    pub fn tlog_zero_mass<T: Num>(q: R, om: T, p: R) -> C {
        ((r_0_same_mass(om, p, q) + q * 2. * p) / (r_0_same_mass(om, p, q) - q * 2. * p)).ln()
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
    pub fn j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        let q2 = q * q;
        let en = energy(q, m);
        bose_distribution_zero_chempot(en, beta) * q2 / (en * 2. * PI2)
    }

    #[inline(always)]
    pub fn j_0_i<T: Num>(q: R, beta: R) -> T {
        Into::<T>::into(bose_distribution_zero_chempot(Into::<T>::into(q), beta) * q / (2. * PI2))
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
        Into::<T>::into(-bose_distribution_zero_chempot(Into::<T>::into(q), beta) * q3 / (2. * PI2))
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
        bose_distribution_zero_chempot(Into::<T>::into(q), beta) * q3 / (6. * PI2)
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
}

#[cfg(test)]
mod tests {
    use crate::low_level::oneloop::thermal;
    use crate::{Num, C, I, R};

    const TOLERANCE: R = 1e-12;

    fn assert_equal<T: Num>(lhs: T, rhs: T) {
        assert!(
            (lhs - rhs).abs() < TOLERANCE,
            "|lhs-rhs| = |({lhs}) - ({rhs})| >= {TOLERANCE:e}"
        );
    }

    fn assert_equal_tol<T: Num>(lhs: T, rhs: T, tol: R) {
        assert!(
            (lhs - rhs).abs() < tol,
            "|lhs-rhs| = |({lhs}) - ({rhs})| >= {tol:e}"
        );
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
                assert_equal_tol(tlog_l(*q, *om, *p, *en, *delta), res[i], 1e-12)
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, en, delta))| {
                assert_equal_tol(oneloop__tlog_l(*q, *om, *p, *en, *delta), res[i], 1e-6)
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
            assert_equal_tol(tlog_l_same_mass(*q, *om, *p, *en), res[i], 1e-6)
        });

        args.iter().enumerate().for_each(|(i, (q, om, p, en))| {
            assert_equal_tol(oneloop__tlog_l_same_mass(*q, *om, *p, *en), res[i], 1e-6)
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

        args.iter().enumerate().for_each(|(i, (q, om, p))| {
            assert_equal_tol(tlog_l_zero_mass(*q, *om, *p), res[i], 1e-6)
        });

        args.iter().enumerate().for_each(|(i, (q, om, p))| {
            assert_equal_tol(oneloop__tlog_l_zero_mass(*q, *om, *p), res[i], 1e-6)
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
    fn test_i_m_m_i_0_m() {
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
                assert_equal(i_m_m_i(*q, *om, *p, 0., *m, *beta), res[i])
            });

        args.iter()
            .enumerate()
            .for_each(|(i, (q, om, p, m, beta))| {
                assert_equal(oneloop__i_m_m_i(*q, *om, *p, 0., *m, *beta), res[i])
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
}
