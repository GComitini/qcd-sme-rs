mod qcd_ym {
    mod gluon {
        use approx::assert_abs_diff_eq;
        use qcd_sme::qcd::{
            thermal::gluon::{
                propagator_l_landau_w_field_config as qcd_propagator_l_landau,
                propagator_t_landau_w_field_config as qcd_propagator_t_landau,
            },
            FieldConfig,
        };
        use qcd_sme::ym::thermal::gluon::{
            propagator_l_landau as ym_propagator_l_landau,
            propagator_t_landau as ym_propagator_t_landau,
        };

        const EPSILON: f64 = 1E-12;

        #[test]
        fn prop_t_landau() {
            let (nc, mg) = (3, 0.656);
            let config = FieldConfig::new(nc, mg, vec![]);

            let test_data = vec![
                (0.0001, 3.5, 4.55, -0.8),
                (0.01, 3.5, 4.55, -0.8),
                (0.01, 0.1, 4.55, -0.8),
                (0.01, 0.1, 1.35, -0.8),
                (0.01, 0.1, 1.35, 1.8),
            ];

            for (om, p, beta, f0) in test_data {
                assert_abs_diff_eq!(
                    qcd_propagator_t_landau(om, p, beta, 0., f0, &config).re,
                    ym_propagator_t_landau(om, p, mg, beta, f0).re,
                    epsilon = EPSILON
                );
            }
        }

        #[test]
        fn prop_l_landau() {
            let (nc, mg) = (3, 0.656);
            let config = FieldConfig::new(nc, mg, vec![]);

            let test_data = vec![
                (0.0001, 3.5, 4.55, -0.8),
                (0.01, 3.5, 4.55, -0.8),
                (0.01, 0.1, 4.55, -0.8),
                (0.01, 0.1, 1.35, -0.8),
                (0.01, 0.1, 1.35, 1.8),
            ];

            for (om, p, beta, f0) in test_data {
                assert_abs_diff_eq!(
                    qcd_propagator_l_landau(om, p, beta, 0., f0, &config).re,
                    ym_propagator_l_landau(om, p, mg, beta, f0).re,
                    epsilon = EPSILON
                );
            }
        }
    }
}

mod thermal_vac {
    mod gluon {
        use approx::assert_abs_diff_eq;
        use qcd_sme::qcd::gluon::propagator_landau_w_field_config as vac_propagator_landau;
        use qcd_sme::qcd::{
            thermal::gluon::{
                propagator_l_landau_w_field_config as propagator_l_landau,
                propagator_l_zero_temp_landau_w_field_config as propagator_l_zero_temp_landau,
                propagator_t_landau_w_field_config as propagator_t_landau,
                propagator_t_zero_temp_landau_w_field_config as propagator_t_zero_temp_landau,
            },
            FieldConfig,
        };

        const EPSILON: f64 = 1E-10;

        #[test]
        fn prop_t_landau() {
            let (nc, mg, beta) = (3, 0.656, 1E6);
            let config = FieldConfig::new(nc, mg, vec![]);

            let test_data = vec![
                (0.0001, 3.5, -0.8),
                (0.01, 3.5, -0.8),
                (0.01, 0.1, -0.8),
                (0.01, 0.1, -0.8),
                (0.01, 0.1, 1.8),
            ];

            let mg2 = mg * mg;

            for (om, p, f0) in test_data {
                let s = (om * om + p * p) / mg2;
                let vac = vac_propagator_landau(s, f0, &config) / mg2;
                assert_abs_diff_eq!(
                    propagator_t_landau(om, p, beta, 0., f0, &config).re,
                    vac,
                    epsilon = EPSILON,
                );
                assert_abs_diff_eq!(
                    propagator_t_zero_temp_landau(om, p, 0., f0, &config).re,
                    vac,
                    epsilon = EPSILON
                );
            }
        }

        #[test]
        fn prop_l_landau() {
            let (nc, mg, beta) = (3, 0.656, 1E6);
            let config = FieldConfig::new(nc, mg, vec![]);

            let test_data = vec![
                (0.0001, 3.5, -0.8),
                (0.01, 3.5, -0.8),
                (0.01, 0.1, -0.8),
                (0.01, 0.1, -0.8),
                (0.01, 0.1, 1.8),
            ];

            let mg2 = mg * mg;

            for (om, p, f0) in test_data {
                let s = (om * om + p * p) / mg2;
                let vac = vac_propagator_landau(s, f0, &config) / mg2;
                assert_abs_diff_eq!(
                    propagator_l_landau(om, p, beta, 0., f0, &config).re,
                    vac,
                    epsilon = EPSILON,
                );
                assert_abs_diff_eq!(
                    propagator_l_zero_temp_landau(om, p, 0., f0, &config).re,
                    vac,
                    epsilon = EPSILON
                );
            }
        }
    }
}
