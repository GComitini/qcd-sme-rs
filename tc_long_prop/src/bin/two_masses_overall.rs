use qcd_sme::qcd::FieldConfig;
use std::fs;
use std::io::{BufWriter, Write};
use tc_long_prop::more_masses::*;

fn main() {
    let basedir = std::path::Path::new(BASEDIR).join("two_masses_overall");
    /* SETTINGS */
    // In units of m
    let (pmin, pmax, dp) = (0.01, 3., 0.01);

    // All but f0 and nc in GeV
    let m = 0.656;
    let f0 = -0.876;
    let nc = 3;
    let om = 1E-5;
    let renpoint = 4.;

    // In units of m
    let temperatures = [0.05, 0.08, 0.12, 0.15, 0.18, 0.21];
    // In GeV
    let quark_masses = [(0.3, 0.4), (0.3, 0.6), (0.8, 1.), (0.3, 1.)];
    let chemical_potentials = [
        0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
    ];

    let critfilepath = basedir.join("tc.out");

    /* DERIVED SETTINGS */
    let prange = pmax - pmin;
    let n = (prange / dp) as usize;
    let momenta: Vec<R> = (0..=n).map(|k| pmin + dp * (k as R)).collect();

    let renpoint2 = renpoint * renpoint;
    let m2 = m * m;
    let renf = renpoint2 / m2;

    fs::create_dir_all(&basedir)
        .unwrap_or_else(|_| panic!("Could not create directory {}", basedir.to_string_lossy()));
    let mut critfile = BufWriter::new(
        fs::File::create(&critfilepath)
            .unwrap_or_else(|_| panic!("Could not create {}", critfilepath.to_string_lossy())),
    );

    // nf = 0
    let mut legends = vec![String::from("$T/m=0.00$")];
    temperatures
        .iter()
        .for_each(|t| legends.push(format!("$T/m={t:.3}$")));

    /* LONGITUDINAL PROPAGATOR AS A FUNCTION OF MOMENTUM,
    MULTIPLE TEMPERATURES */
    println!("*** Plotting mq1 = N/A GeV, mq2 = N/A, mu = N/A GeV, nf = 0 ***");
    let mut plot = Plot2D::new();
    plot.set_domain(momenta.clone());
    plot.set_ylim((0., 5.));
    plot.set_xlabel("$p/m$");
    plot.set_ylabel("$m^{2}\\Delta_{L}(p)$");
    plot.set_title("$M_{q}^{(1)}=$N/A GeV, $M_{q}^{(2)}=$N/A GeV, $n_{F}=0$, $\\mu=$N/A GeV");

    // T=0
    let z = ym_propagator_l_zero_temp_landau((om * om + renpoint * renpoint) / m2, f0) * renf;
    println!("Computed renormalization factor: z = {z}");
    plot.insert_image(
        momenta
            .iter()
            .map(|&p| {
                println!(
                    "Computing p/m = {p:.3} @ (T/m, mq1, mq2, mu) = (0.00, N/A GeV, N/A GeV, N/A GeV), nf = 0..."
                );
                ym_propagator_l_zero_temp_landau(om * om / m2 + p * p, f0) / z
            })
            .collect(),
    );

    // T>0
    for t in temperatures {
        let beta = 1. / (t * m);
        let z = ym_propagator_l_landau(om, renpoint, m, beta, f0).re * renf;
        println!("Computed renormalization factor: z = {z}");
        plot.insert_image(
                        momenta
                            .iter()
                            .map(|&p| {
                                println!(
                                "Computing p/m = {p:.3} @ (T/m, mq1, mq2, mu) = ({t:.3}, N/A GeV, N/A GeV, N/A GeV), nf = 0..."
                            );
                                ym_propagator_l_landau(om, p * m, m,  beta,  f0).re / z
                            })
                            .collect(),
                    );
    }

    plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
    plot.set_path(&format!(
        "{}/nf_0_mq1_NA_mq2_NA_mu_NA.png",
        basedir.to_string_lossy()
    ));
    plot.savefig().expect("Could not save figure");

    /* LONGITUDINAL PROPAGATOR AT SMALL MOMENTUM AS A FUNCTION OF TEMPERATURE */
    let pbase = pmin;
    let (tmin, tmax, dt) = (0.01, 0.2, 0.005);
    let trange = tmax - tmin;
    let n = (trange / dt) as usize;
    let tdomain: Vec<R> = (0..=n).map(|k| tmin + dt * (k as R)).collect();

    println!("*** Plotting mq1 = N/A GeV, mq2 = N/A GeV, nf = 0 ***");
    let mut plot = Plot2D::new();
    let mut domain = vec![0.];
    domain.extend(tdomain.clone());
    plot.set_domain(domain.clone());
    plot.set_ylim((0., 5.));
    plot.set_xlabel("$T/m$");
    plot.set_ylabel(&format!("$m^{{2}}\\Delta_{{L}}(p={pbase}\\,m)$"));
    plot.set_title("$M_{q}^{(1)}=$ N/A GeV, $M_{q}^{(2)}=$ N/A GeV, $n_{F}=0$");

    // T=0
    let z = ym_propagator_l_zero_temp_landau((om * om + renpoint * renpoint) / m2, f0) * renf;
    println!("Computed renormalization factor: z = {z}");

    println!("Computing p/m = {pbase:.3} @ (T/m, mq1, mq2, mu) = (0.00, N/A GeV, N/A GeV,N/A GeV), nf = 0...");
    let mut data = vec![ym_propagator_l_zero_temp_landau(om * om / m2 + pbase * pbase, f0) / z];

    // T>0
    tdomain.iter().for_each(|t| {
        let beta = 1. / (t * m);
        let z = ym_propagator_l_landau(om, renpoint, m, beta, f0).re * renf;
        println!("Computed renormalization factor: z = {z}");
        println!(
            "Computing p/m = {pbase:.3} @ (T/m, mq1, mq2, mu) = ({t:.3}, N/A GeV, N/A GeV, N/A GeV), nf = 0..."
        );
        data.push(ym_propagator_l_landau(om, pbase * m, m, beta, f0).re / z);
    });

    let tc = find_critical_temperature(&domain, &data);
    writeln!(critfile, "0\t0.000\t0.000\t0.000\t{tc:.3}")
        .unwrap_or_else(|_| panic!("Could not write to {}", critfilepath.to_string_lossy()));

    plot.insert_image(data);

    plot.set_path(&format!(
        "{}/ir_nf_0_mq1_NA_mq2_NA.png",
        basedir.to_string_lossy()
    ));
    plot.savefig().expect("Could not save figure");

    // nf > 0
    for nf in [1, 4] {
        let mut legends = vec![String::from("$T/m=0.00$")];
        temperatures
            .iter()
            .for_each(|t| legends.push(format!("$T/m={t:.3}$")));

        /* LONGITUDINAL PROPAGATOR AS A FUNCTION OF MOMENTUM
        AT FIXED MQ, MU AND NF, MULTIPLE TEMPERATURES */
        for mq in quark_masses {
            let (mq1, mq2) = (mq.0, mq.1);
            let config = FieldConfig::new(nc, m, vec![(nf, mq1), (nf, mq2)]);
            for mu in chemical_potentials {
                println!("*** Plotting mq = {mq1:.3} GeV, mq = {mq2:.3} GeV, mu = {mu:.3} GeV, nf = {nf} ***");
                let mut plot = Plot2D::new();
                plot.set_domain(momenta.clone());
                plot.set_ylim((0., 5.));
                plot.set_xlabel("$p/m$");
                plot.set_ylabel("$m^{2}\\Delta_{L}(p)$");
                plot.set_title(&format!(
                    "$M_{{q}}^{{(1)}}={mq1:.3}$ GeV, $M_{{q}}^{{(2)}}={mq2:.3}$ GeV, $n_{{F}}={nf}$, $\\mu={mu:.3}$ GeV"
                ));

                // T=0
                let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, &config).re * renf;
                println!("Computed renormalization factor: z = {z}");
                plot.insert_image(
                    momenta
                        .iter()
                        .map(|&p| {
                            println!(
                            "Computing p/m = {p:.3} @ (T/m, mq1, mq2, mu) = (0.00, {mq1:.3} GeV, {mq2:.3} GeV, {mu:.3} GeV), nf = {nf}..."
                        );
                            propagator_l_zero_temp_landau(om, p * m,  mu, f0, &config).re / z
                        })
                        .collect(),
                );

                // T>0
                for t in temperatures {
                    let beta = 1. / (t * m);
                    let z = propagator_l_landau(om, renpoint, beta, mu, f0, &config).re * renf;
                    println!("Computed renormalization factor: z = {z}");
                    plot.insert_image(
                        momenta
                            .par_iter()
                            .map(|&p| {
                                println!(
                                "Computing p/m = {p:.3} @ (T/m, mq1, mq2, mu) = ({t:.3}, {mq1:.3} GeV, {mq2:.3} GeV, {mu:.3} GeV), nf = {nf}..."
                            );
                                propagator_l_landau(om, p * m, beta, mu, f0, &config).re / z
                            })
                            .collect(),
                    );
                }

                plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
                plot.set_path(&format!(
                    "{}/nf_{nf}_mq1_{mq1:.3}_mq2_{mq2:.3}_mu_{mu:.3}.png",
                    basedir.to_string_lossy()
                ));
                plot.savefig().expect("Could not save figure");
            }
        }

        /* LONGITUDINAL PROPAGATOR AT SMALL MOMENTUM AS A FUNCTION OF TEMPERATURE
        AT FIXED MQ AND NF, MULTIPLE CHEMICAL POTENTIALS */
        let pbase = pmin;
        let (tmin, tmax, dt) = (0.01, 0.2, 0.005);
        let trange = tmax - tmin;
        let n = (trange / dt) as usize;
        let tdomain: Vec<R> = (0..=n).map(|k| tmin + dt * (k as R)).collect();
        let legends: Vec<String> = chemical_potentials
            .iter()
            .map(|mu| format!("$\\mu={mu:.3}$ GeV"))
            .collect();

        for mq in quark_masses {
            let (mq1, mq2) = (mq.0, mq.1);
            let config = FieldConfig::new(nc, m, vec![(nf, mq1), (nf, mq2)]);
            println!("*** Plotting mq1 = {mq1:.3} GeV, mq2 = {mq2:.3} GeV, nf = {nf} ***");
            let mut plot = Plot2D::new();
            let mut domain = vec![0.];
            domain.extend(tdomain.clone());
            plot.set_domain(domain.clone());
            plot.set_ylim((0., 5.));
            plot.set_xlabel("$T/m$");
            plot.set_ylabel(&format!("$m^{{2}}\\Delta_{{L}}(p={pbase}\\,m)$"));
            plot.set_title(&format!(
                "$M_{{q}}^{{(1)}}={mq1:.3}$ GeV, $M_{{q}}^{{(2)}}={mq2:.3}$ GeV, $n_{{F}}={nf}$"
            ));

            // Also collect critical temperatures for later plotting
            let mut tcs = vec![];

            for mu in chemical_potentials {
                // T=0
                let z = propagator_l_zero_temp_landau(om, renpoint, mu, f0, &config).re * renf;
                println!("Computed renormalization factor: z = {z}");

                println!("Computing p/m = {pbase:.3} @ (T/m, mq1, mq2, mu) = (0.00, {mq1:.3} GeV, {mq2:.3} GeV, {mu:.3} GeV), nf = {nf}...");
                let mut data =
                    vec![propagator_l_zero_temp_landau(om, pbase * m, mu, f0, &config).re / z];

                // T>0
                tdomain.iter().for_each(|t| {
                    let beta = 1./(t*m);
                    let z = propagator_l_landau(om, renpoint, beta, mu, f0, &config).re * renf;
                    println!("Computed renormalization factor: z = {z}");
                    println!(
                        "Computing p/m = {pbase:.3} @ (T/m, mq1, mq2, mu) = ({t:.3}, {mq1:.3} GeV, {mq2:.3} GeV, {mu:.3} GeV), nf = {nf}...");
                    data.push(propagator_l_landau(om, pbase * m, beta, mu, f0, &config).re / z);
                });

                let tc = find_critical_temperature(&domain, &data);
                tcs.push(tc);
                writeln!(critfile, "{nf}\t{mq1:.3}\t{mq2:.3}\t{mu:.3}\t{tc:.3}").unwrap_or_else(
                    |_| panic!("Could not write to {}", critfilepath.to_string_lossy()),
                );

                plot.insert_image(data);
            }

            plot.set_legend(legends.iter().map(|s| s.as_str()).collect());
            plot.set_path(&format!(
                "{}/ir_nf_{nf}_mq1_{mq1:.3}_mq2_{mq2:.3}.png",
                basedir.to_string_lossy()
            ));
            plot.savefig().expect("Could not save figure");

            /* CRITICAL TEMPERATURE AS A FUNCTION OF THE CHEMICAL POTENTIAL
            AT FIXED MQ AND NF */
            let mut plot = Plot2D::new();
            plot.set_domain(chemical_potentials.to_vec());
            plot.set_ylim((0., 0.22));
            plot.set_xlabel("$\\mu$ [GeV]");
            plot.set_ylabel("$T_{c}/m$");
            plot.set_title(&format!(
                "$M_{{q}}^{{(1)}}={mq1:.3}$ GeV, $M_{{q}}^{{(2)}}={mq2:.3}$ GeV, $n_{{F}}={nf}$"
            ));
            plot.insert_image(tcs);
            plot.set_path(&format!(
                "{}/tc_nf_{nf}_mq1_{mq1:.3}_mq2_{mq2:.3}.png",
                basedir.to_string_lossy()
            ));
            plot.savefig().expect("Could not save figure");
        }
    }
}
