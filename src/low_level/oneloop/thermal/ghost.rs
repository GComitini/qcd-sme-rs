// For use in other languages, e.g. C/C++/Python
//
// - Re-export at crate::ffi, since symbols need to be unmangled anyway and the
//   namespace will not be preserved.
//
// - Here too each function calls its inline analogue, but the objective is to
//   switch from generic to concrete argument types so that the functions can
//   be compiled into a C dynamic library. To do so, we need to double their
//   number (one function for real arguments, another for complex arguments).
pub(crate) mod ffi {}

// For internal use only
//
// - Here we define the building blocks for the other functions. This module
//   serves two purposes: to hold inlined functions and to provide a single
//   source of truth for the actual mathematical expressions
pub(crate) mod inlines {
    // See if it's better not to inline these
    use crate::low_level::oneloop::thermal::*;
    use crate::{Num, C, R};

    #[inline(always)]
    pub fn thermal_self_energy_landau_i<T: Num>(q: R, om: T, p: R, m: R, beta: R) -> C {
        let s = om * om + p * p;
        let m2 = m * m;
        let a = s + m2;

        let t1 = -a * s / (2. * m2) * i_m_0_i(q, om, p, m, beta);
        let t2 = s * s / (2. * m2) * i_0_0_i(q, om, p, beta);
        let t3 = (s * 2. - m2) / (4. * m2) * (j_m_i(q, m, beta) - j_0_i(q, beta));
        let t4 = a * a / 4. * d_i_m_0_i(q, om, p, m, beta);
        let t5 = -(s - m2) / 4. * d_j_m_i(q, m, beta);

        t5 + (t3 + t1) + t2 + t4
    }
}
