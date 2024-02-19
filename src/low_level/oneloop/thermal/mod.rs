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

    pub fn r_0<T1: Num, T2: Num>(om: T1, p: R, en: T2, delta: R) -> C {
        inlines::r_0(om, p, en, delta)
    }

    pub fn tlog<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2, delta: R) -> C {
        inlines::tlog(q, om, p, en, delta)
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
}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    use crate::common::{energy, thermal::bose_distribution_zero_chempot};
    use crate::{Num, C, I, R};
    use std::f64::consts::PI;

    const PI2: R = PI * PI;

    #[inline(always)]
    pub fn r_0<T1: Num, T2: Num>(om: T1, p: R, en: T2, delta: R) -> C {
        om * om + om * (en * 2. * I) + p * p + delta
    }

    #[inline(always)]
    pub fn tlog<T1: Num, T2: Num>(q: R, om: T1, p: R, en: T2, delta: R) -> C {
        ((r_0(om, p, en, delta) + q * 2. * p) / (r_0(om, p, en, delta) - q * 2. * p)).ln()
    }

    #[inline(always)]
    pub fn j_m_i<T: Num>(q: R, m: T, beta: R) -> T {
        let q2 = q * q;
        let en = energy(q, m);
        bose_distribution_zero_chempot(en, beta) * q2 / (en * 2. * PI2)
    }

    #[inline(always)]
    pub fn j_0_i<T: Num>(q: R, beta: R) -> T {
        (bose_distribution_zero_chempot(Into::<T>::into(q), beta) * q / (2. * PI2)).into()
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
        (-bose_distribution_zero_chempot(Into::<T>::into(q), beta) * q3 / (2. * PI2)).into()
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

    #[test]
    fn test_tlog() {
        use thermal::{ffi::oneloop__tlog, tlog};

        let args: [(R, R, R, R, R); 6] = [
            (0.03, 0., 2.12, 4., 2.),
            (0.03, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 2.12, 4., 2.),
            (0.047, 0.18, 1.27, 4., 2.),
            (0.047, 0.18, 1.27, 6.9, 2.),
            (0.047, 0.18, 1.27, 6.9, 3.21),
        ];

        let res: [C; 6] = [
            0.039177220079546576 + 0. * I,
            0.03717215443439166 - 0.008203228369481758 * I,
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
}
