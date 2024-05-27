use crate::types::{NCTYPE, NFTYPE};
use crate::R;

/// Field configuration for gluons and quarks.
pub struct FieldConfig {
    /// Number of colors.
    pub nc: NCTYPE,
    /// Gluon mass parameter.
    pub gluon: R,
    /// Quark mass configuration: number of quarks (`quarks[i].0`) with a given
    /// mass parameter (`quarks[i].1`).
    pub quarks: Vec<(NFTYPE, R)>,
    // (nf, mq, nf/nc, mq/m)
    pub(crate) quarks_internal: Vec<QuarkInternal>,
}

impl FieldConfig {
    pub fn new(nc: NCTYPE, gluon: R, quarks: Vec<(NFTYPE, R)>) -> Self {
        let quarks_internal = quarks
            .iter()
            .map(|&(nf, mq)| QuarkInternal::new(nc, gluon, nf, mq))
            .collect();
        Self {
            nc,
            gluon,
            quarks,
            quarks_internal,
        }
    }
}

pub struct QuarkInternal {
    pub nf: NFTYPE,
    pub mq: R,
    pub nf_div_nc: R,
    pub mq_div_m_2: R,
}

impl QuarkInternal {
    pub(crate) fn new(nc: NCTYPE, m: R, nf: NFTYPE, mq: R) -> Self {
        let nf_div_nc = (nf as R) / (nc as R);
        let mq_div_m = mq / m;
        let mq_div_m_2 = mq_div_m * mq_div_m;
        Self {
            nf,
            mq,
            nf_div_nc,
            mq_div_m_2,
        }
    }
}
