use qcd_sme::{
    qcd::FieldConfig,
    types::{NCTYPE, R},
};

const NC: NCTYPE = 3;

fn main() {
    [
        (0.139, 0.7360461258894254, -0.38326508306172924),
        (0.154, 0.7492957152865027, -0.31556760023282565),
        (0.174, 0.7217297548149253, -0.25830998457395093),
        (0.199, 0.7196745401373572, -0.24971877005975093),
        (0.233, 0.7375872973499797, -0.21776224853123094),
        (0.278, 0.7855533106315216, -0.1736425033681984),
    ]
    .iter()
    .for_each(|&(t, mg, f00)| {
        let fieldconfig = FieldConfig::new(NC, mg, vec![(2, 0.2)]);
        let f0 = f00
            - fieldconfig
                .quarks
                .iter()
                .map(|(nf, mq)| (*nf as R) * (fieldconfig.gluon / mq).ln())
                .sum::<R>()
                * 4.
                / (9. * (fieldconfig.nc as R));
        println!("{:.0} & {:.1} & {f0:.4}\\\\", t * 1000., mg * 1000.);
    });
}
