use qcd_sme::{
    qcd::FieldConfig,
    types::{NCTYPE, R},
};

const NC: NCTYPE = 3;

fn main() {
    [
        (0.139, 0.7517777951910817, -0.5062071238487814),
        (0.154, 0.7639062870106536, -0.44084616449186326),
        (0.174, 0.7345150400347596, -0.3843769966635342),
        (0.199, 0.7306060280562273, -0.38161604237819363),
        (0.233, 0.746307345288752, -0.3585815488317115),
        (0.278, 0.7924693589714722, -0.32681035700549366),
    ]
    .iter()
    .for_each(|&(t, mg, f00)| {
        let fieldconfig = FieldConfig::new(NC, mg, vec![(2, 0.4)]);
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
