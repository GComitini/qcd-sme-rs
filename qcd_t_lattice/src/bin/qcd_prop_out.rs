use qcd_sme::consts::{set_default_max_iter_integral, set_default_tol_integral};
use qcd_t_lattice::{gluon::*, FieldConfig, OMEPS, R};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

const PMIN: R = 0.01;
const PMAX: R = 4.;
const NP: usize = 400;
const DP: R = (PMAX - PMIN) / (NP as R);

fn main() {
    set_default_max_iter_integral(30);
    set_default_tol_integral(1E-7);

    let comp = std::env::args()
        .nth(1)
        .expect("You need to pass the component as argument!");

    match comp.as_str() {
        "l" | "L" | "t" | "T" => {}
        _ => panic!("The component argument can only be t, T, l or L!"),
    }

    let t = std::env::args()
        .nth(2)
        .expect("You need to pass the temperature in GeV as argument!")
        .parse::<R>()
        .expect("The temperature argument must be a float!");

    let mq = std::env::args()
        .nth(3)
        .expect("You need to pass the quark mass in GeV as argument!")
        .replace(",", "")
        .parse::<R>()
        .expect("The quark mass argument must be a float!");

    let mg = std::env::args()
        .nth(4)
        .expect("You need to pass the gluon mass in GeV as argument!")
        .replace(",", "")
        .parse::<R>()
        .expect("The gluon mass argument must be a float!");

    let z = std::env::args()
        .nth(5)
        .expect("You need to pass z as argument!")
        .replace(",", "")
        .parse::<R>()
        .expect("The z argument must be a float!");

    let f0 = std::env::args()
        .nth(6)
        .expect("You need to pass f0 as argument!")
        .parse::<R>()
        .expect("The f0 argument must be a float!");

    let ff = &|p| {
        let fieldconfig = FieldConfig::new(3, mg, vec![(2, mq)]);
        match comp.as_str() {
            "l" | "L" => z * propagator_l(OMEPS, p, 1. / t, f0, &fieldconfig).re,
            "t" | "T" => z * propagator_t(OMEPS, p, 1. / t, f0, &fieldconfig).re,
            _ => unreachable!(),
        }
    };

    let ps: Vec<R> = (0..=NP).map(|i| PMIN + DP * (i as R)).collect();

    let prop: Vec<(R, R)> = ps.par_iter().map(|p| (*p, ff(*p))).collect();

    let mut fout = BufWriter::new(
        File::create(format!(
            "{}{:.0}mq{}mg{}.txt",
            comp.to_uppercase(),
            (t * 1000.) as u64,
            (mq * 1000.) as u64,
            (mg * 1000.) as u64
        ))
        .unwrap(),
    );

    prop.iter()
        .for_each(|(p, pp)| writeln!(fout, "{p} {pp}").unwrap());
}
