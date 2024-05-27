mod qcd {
    mod gluon {
        use qcd_sme::consts::{set_number_of_colors, set_number_of_fermions};
        use qcd_sme::qcd::FieldConfig;
        use qcd_sme::{Num, C, R};
        use serial_test::serial;

        const TOLERANCE: R = 1e-12;

        const REAL_TEST_VAL: [R; 4] = [1., 1.834, 2.5, 8.18];
        const COMPLEX_TEST_VAL: [C; 5] = [
            C { re: 1., im: -0.5 },
            C { re: 1., im: 0.5 },
            C { re: 1.347, im: 0. },
            C {
                re: 2.56,
                im: 0.732,
            },
            C {
                re: 7.28,
                im: 5.166,
            },
        ];

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
        #[serial]
        fn test_prop_landau_w_field_config() {
            use qcd_sme::qcd::gluon::dressing_inv_landau;
            use qcd_sme::qcd::gluon::{
                dressing_landau_w_field_config, propagator_landau_w_field_config,
            };
            use qcd_sme::ym::gluon::dressing_inv_landau as ym_dressing_inv_landau;

            let mg = 0.656;

            let ncs = [2, 3];

            for nc in ncs {
                let config = FieldConfig::new(nc, mg, vec![]);

                let real_results: Vec<R> = REAL_TEST_VAL
                    .iter()
                    .map(|&s| ym_dressing_inv_landau(s, -0.876))
                    .collect();
                let complex_results: Vec<C> = COMPLEX_TEST_VAL
                    .iter()
                    .map(|&s| ym_dressing_inv_landau(s, -0.876))
                    .collect();

                REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                    assert_equal(
                        1. / dressing_landau_w_field_config(s, -0.876, &config),
                        real_results[i],
                    );
                    assert_equal(
                        1. / (s * propagator_landau_w_field_config(s, -0.876, &config)),
                        real_results[i],
                    )
                });

                COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                    assert_equal(
                        1. / dressing_landau_w_field_config(s, -0.876, &config),
                        complex_results[i],
                    );
                    assert_equal(
                        1. / (s * propagator_landau_w_field_config(s, -0.876, &config)),
                        complex_results[i],
                    )
                });
            }

            for nc in ncs {
                set_number_of_colors(nc);
                for nq in [1, 2] {
                    set_number_of_fermions(nq);
                    for mq in [0.1, 0.5] {
                        let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);

                        let real_results: Vec<R> = REAL_TEST_VAL
                            .iter()
                            .map(|&s| dressing_inv_landau(s, mq / mg, -0.876))
                            .collect();
                        let complex_results: Vec<C> = COMPLEX_TEST_VAL
                            .iter()
                            .map(|&s| dressing_inv_landau(s, mq / mg, -0.876))
                            .collect();

                        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                            assert_equal(
                                1. / dressing_landau_w_field_config(s, -0.876, &config),
                                real_results[i],
                            );
                            assert_equal(
                                1. / (s * propagator_landau_w_field_config(s, -0.876, &config)),
                                real_results[i],
                            )
                        });

                        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                            assert_equal(
                                1. / dressing_landau_w_field_config(s, -0.876, &config),
                                complex_results[i],
                            );
                            assert_equal(
                                1. / (s * propagator_landau_w_field_config(s, -0.876, &config)),
                                complex_results[i],
                            )
                        });
                    }
                }
            }

            set_number_of_colors(3);
            set_number_of_fermions(1);
        }

        #[test]
        #[serial]
        fn test_prop_w_field_config() {
            use qcd_sme::qcd::gluon::dressing_inv;
            use qcd_sme::qcd::gluon::{dressing_w_field_config, propagator_w_field_config};
            use qcd_sme::ym::gluon::dressing_inv as ym_dressing_inv;

            let mg = 0.656;

            let ncs = [2, 3];
            let xis = [0., 1.];

            for xi in xis {
                for nc in ncs {
                    let config = FieldConfig::new(nc, mg, vec![]);

                    let real_results: Vec<R> = REAL_TEST_VAL
                        .iter()
                        .map(|&s| ym_dressing_inv(s, -0.876, xi))
                        .collect();
                    let complex_results: Vec<C> = COMPLEX_TEST_VAL
                        .iter()
                        .map(|&s| ym_dressing_inv(s, -0.876, xi))
                        .collect();

                    REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                        assert_equal(
                            1. / dressing_w_field_config(s, -0.876, xi, &config),
                            real_results[i],
                        );
                        assert_equal(
                            1. / (s * propagator_w_field_config(s, -0.876, xi, &config)),
                            real_results[i],
                        )
                    });

                    COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                        assert_equal(
                            1. / dressing_w_field_config(s, -0.876, xi, &config),
                            complex_results[i],
                        );
                        assert_equal(
                            1. / (s * propagator_w_field_config(s, -0.876, xi, &config)),
                            complex_results[i],
                        )
                    });
                }
            }

            for xi in xis {
                for nc in ncs {
                    set_number_of_colors(nc);
                    for nq in [1, 2] {
                        set_number_of_fermions(nq);
                        for mq in [0.1, 0.5] {
                            let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);

                            let real_results: Vec<R> = REAL_TEST_VAL
                                .iter()
                                .map(|&s| dressing_inv(s, mq / mg, -0.876, xi))
                                .collect();
                            let complex_results: Vec<C> = COMPLEX_TEST_VAL
                                .iter()
                                .map(|&s| dressing_inv(s, mq / mg, -0.876, xi))
                                .collect();

                            REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                                assert_equal(
                                    1. / dressing_w_field_config(s, -0.876, xi, &config),
                                    real_results[i],
                                );
                                assert_equal(
                                    1. / (s * propagator_w_field_config(s, -0.876, xi, &config)),
                                    real_results[i],
                                )
                            });

                            COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                                assert_equal(
                                    1. / dressing_w_field_config(s, -0.876, xi, &config),
                                    complex_results[i],
                                );
                                assert_equal(
                                    1. / (s * propagator_w_field_config(s, -0.876, xi, &config)),
                                    complex_results[i],
                                )
                            });
                        }
                    }
                }
            }
            set_number_of_colors(3);
            set_number_of_fermions(1);
        }

        #[test]
        #[serial]
        fn test_prop_crossed_landau_w_field_config() {
            use qcd_sme::qcd::gluon::dressing_crossed_inv_landau;
            use qcd_sme::qcd::gluon::{
                dressing_crossed_landau_w_field_config, propagator_crossed_landau_w_field_config,
            };
            use qcd_sme::ym::gluon::dressing_inv_landau as ym_dressing_inv_landau;

            let mg = 0.656;

            let ncs = [2, 3];

            for nc in ncs {
                let config = FieldConfig::new(nc, mg, vec![]);

                let real_results: Vec<R> = REAL_TEST_VAL
                    .iter()
                    .map(|&s| ym_dressing_inv_landau(s, -0.876))
                    .collect();
                let complex_results: Vec<C> = COMPLEX_TEST_VAL
                    .iter()
                    .map(|&s| ym_dressing_inv_landau(s, -0.876))
                    .collect();

                REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                    assert_equal(
                        1. / dressing_crossed_landau_w_field_config(s, -0.876, &config),
                        real_results[i],
                    );
                    assert_equal(
                        1. / (s * propagator_crossed_landau_w_field_config(s, -0.876, &config)),
                        real_results[i],
                    )
                });

                COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                    assert_equal(
                        1. / dressing_crossed_landau_w_field_config(s, -0.876, &config),
                        complex_results[i],
                    );
                    assert_equal(
                        1. / (s * propagator_crossed_landau_w_field_config(s, -0.876, &config)),
                        complex_results[i],
                    )
                });
            }

            for nc in ncs {
                set_number_of_colors(nc);
                for nq in [1, 2] {
                    set_number_of_fermions(nq);
                    for mq in [0.1, 0.5] {
                        let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);

                        let real_results: Vec<R> = REAL_TEST_VAL
                            .iter()
                            .map(|&s| dressing_crossed_inv_landau(s, mq / mg, -0.876))
                            .collect();
                        let complex_results: Vec<C> = COMPLEX_TEST_VAL
                            .iter()
                            .map(|&s| dressing_crossed_inv_landau(s, mq / mg, -0.876))
                            .collect();

                        REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                            assert_equal(
                                1. / dressing_crossed_landau_w_field_config(s, -0.876, &config),
                                real_results[i],
                            );
                            assert_equal(
                                1. / (s * propagator_crossed_landau_w_field_config(
                                    s, -0.876, &config,
                                )),
                                real_results[i],
                            )
                        });

                        COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                            assert_equal(
                                1. / dressing_crossed_landau_w_field_config(s, -0.876, &config),
                                complex_results[i],
                            );
                            assert_equal(
                                1. / (s * propagator_crossed_landau_w_field_config(
                                    s, -0.876, &config,
                                )),
                                complex_results[i],
                            )
                        });
                    }
                }
            }
            set_number_of_colors(3);
            set_number_of_fermions(1);
        }

        #[test]
        #[serial]
        fn test_prop_crossed_w_field_config() {
            use qcd_sme::qcd::gluon::dressing_crossed_inv;
            use qcd_sme::qcd::gluon::{
                dressing_crossed_w_field_config, propagator_crossed_w_field_config,
            };
            use qcd_sme::ym::gluon::dressing_inv as ym_dressing_inv;

            let mg = 0.656;

            let ncs = [2, 3];
            let xis = [0., 1.];

            for xi in xis {
                for nc in ncs {
                    let config = FieldConfig::new(nc, mg, vec![]);

                    let real_results: Vec<R> = REAL_TEST_VAL
                        .iter()
                        .map(|&s| ym_dressing_inv(s, -0.876, xi))
                        .collect();
                    let complex_results: Vec<C> = COMPLEX_TEST_VAL
                        .iter()
                        .map(|&s| ym_dressing_inv(s, -0.876, xi))
                        .collect();

                    REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                        assert_equal(
                            1. / dressing_crossed_w_field_config(s, -0.876, xi, &config),
                            real_results[i],
                        );
                        assert_equal(
                            1. / (s * propagator_crossed_w_field_config(s, -0.876, xi, &config)),
                            real_results[i],
                        )
                    });

                    COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                        assert_equal(
                            1. / dressing_crossed_w_field_config(s, -0.876, xi, &config),
                            complex_results[i],
                        );
                        assert_equal(
                            1. / (s * propagator_crossed_w_field_config(s, -0.876, xi, &config)),
                            complex_results[i],
                        )
                    });
                }
            }

            for xi in xis {
                for nc in ncs {
                    set_number_of_colors(nc);
                    for nq in [1, 2] {
                        set_number_of_fermions(nq);
                        for mq in [0.1, 0.5] {
                            let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);

                            let real_results: Vec<R> = REAL_TEST_VAL
                                .iter()
                                .map(|&s| dressing_crossed_inv(s, mq / mg, -0.876, xi))
                                .collect();
                            let complex_results: Vec<C> = COMPLEX_TEST_VAL
                                .iter()
                                .map(|&s| dressing_crossed_inv(s, mq / mg, -0.876, xi))
                                .collect();

                            REAL_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                                assert_equal(
                                    1. / dressing_crossed_w_field_config(s, -0.876, xi, &config),
                                    real_results[i],
                                );
                                assert_equal(
                                    1. / (s * propagator_crossed_w_field_config(
                                        s, -0.876, xi, &config,
                                    )),
                                    real_results[i],
                                )
                            });

                            COMPLEX_TEST_VAL.iter().enumerate().for_each(|(i, &s)| {
                                assert_equal(
                                    1. / dressing_crossed_w_field_config(s, -0.876, xi, &config),
                                    complex_results[i],
                                );
                                assert_equal(
                                    1. / (s * propagator_crossed_w_field_config(
                                        s, -0.876, xi, &config,
                                    )),
                                    complex_results[i],
                                )
                            });
                        }
                    }
                }
            }
            set_number_of_colors(3);
            set_number_of_fermions(1);
        }
    }

    mod thermal {
        mod gluon {
            use qcd_sme::consts::{set_number_of_colors, set_number_of_fermions};
            use qcd_sme::{Num, R};
            use serial_test::serial;

            const TOLERANCE: R = 1e-12;

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
            #[serial]
            fn test_prop_l_w_field_config() {
                use qcd_sme::qcd::thermal::gluon::{
                    propagator_l_landau, propagator_l_landau_w_field_config,
                };
                use qcd_sme::qcd::FieldConfig;
                use qcd_sme::ym::thermal::gluon::propagator_l_landau as ym_propagator_l_landau;

                let mg = 0.656;
                let f0 = -0.876;

                let ncs = [2, 3];
                let nqs = [1, 2];
                let mqs = [0.1, 0.5];
                let ts = [0.05, 0.4];
                let mus = [0., 0.7];
                let oms = [0.37, 1.13];
                let ps = [0.41, 0.94];

                for nc in ncs {
                    set_number_of_colors(nc);
                    let config = FieldConfig::new(nc, mg, vec![]);
                    for t in ts {
                        let beta = 1. / t;
                        for om in oms {
                            for p in ps {
                                assert_equal(
                                    ym_propagator_l_landau(om, p, mg, beta, f0),
                                    propagator_l_landau_w_field_config(
                                        om, p, beta, 0., f0, &config,
                                    ),
                                )
                            }
                        }
                    }
                }

                for nc in ncs {
                    set_number_of_colors(nc);
                    for nq in nqs {
                        set_number_of_fermions(nq);
                        for mq in mqs {
                            let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);
                            for t in ts {
                                let beta = 1. / t;
                                for mu in mus {
                                    for om in oms {
                                        for p in ps {
                                            assert_equal(
                                                propagator_l_landau(om, p, mg, mq, beta, mu, f0),
                                                propagator_l_landau_w_field_config(
                                                    om, p, beta, mu, f0, &config,
                                                ),
                                            )
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                set_number_of_colors(3);
                set_number_of_fermions(1);
            }

            #[test]
            #[serial]
            fn test_prop_t_w_field_config() {
                use qcd_sme::qcd::thermal::gluon::{
                    propagator_t_landau, propagator_t_landau_w_field_config,
                };
                use qcd_sme::qcd::FieldConfig;
                use qcd_sme::ym::thermal::gluon::propagator_t_landau as ym_propagator_t_landau;

                let mg = 0.656;
                let f0 = -0.876;

                let ncs = [2, 3];
                let nqs = [1, 2];
                let mqs = [0.1, 0.5];
                let ts = [0.05, 0.4];
                let mus = [0., 0.7];
                let oms = [0.37, 1.13];
                let ps = [0.41, 0.94];

                for nc in ncs {
                    set_number_of_colors(nc);
                    let config = FieldConfig::new(nc, mg, vec![]);
                    for t in ts {
                        let beta = 1. / t;
                        for om in oms {
                            for p in ps {
                                assert_equal(
                                    ym_propagator_t_landau(om, p, mg, beta, f0),
                                    propagator_t_landau_w_field_config(
                                        om, p, beta, 0., f0, &config,
                                    ),
                                )
                            }
                        }
                    }
                }

                for nc in ncs {
                    set_number_of_colors(nc);
                    for nq in nqs {
                        set_number_of_fermions(nq);
                        for mq in mqs {
                            let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);
                            for t in ts {
                                let beta = 1. / t;
                                for mu in mus {
                                    for om in oms {
                                        for p in ps {
                                            assert_equal(
                                                propagator_t_landau(om, p, mg, mq, beta, mu, f0),
                                                propagator_t_landau_w_field_config(
                                                    om, p, beta, mu, f0, &config,
                                                ),
                                            )
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                set_number_of_colors(3);
                set_number_of_fermions(1);
            }

            #[test]
            #[serial]
            fn test_prop_l_zero_temp_w_field_config() {
                use qcd_sme::qcd::thermal::gluon::{
                    propagator_l_zero_temp_landau, propagator_l_zero_temp_landau_w_field_config,
                };
                use qcd_sme::qcd::FieldConfig;
                use qcd_sme::ym::gluon::propagator_landau as ym_propagator_l_zero_temp_landau;

                let mg = 0.656;
                let f0 = -0.876;

                let ncs = [2, 3];
                let nqs = [1, 2];
                let mqs = [0.1, 0.5];
                let mus = [0., 0.7];
                let oms = [0.37, 1.13];
                let ps = [0.41, 0.94];

                let mg2 = mg * mg;

                for nc in ncs {
                    set_number_of_colors(nc);
                    let config = FieldConfig::new(nc, mg, vec![]);
                    for om in oms {
                        for p in ps {
                            let s = (om * om + p * p) / mg2;
                            assert_equal(
                                ym_propagator_l_zero_temp_landau(s, f0) / mg2,
                                propagator_l_zero_temp_landau_w_field_config(
                                    om, p, 0., f0, &config,
                                )
                                .re,
                            )
                        }
                    }
                }

                for nc in ncs {
                    set_number_of_colors(nc);
                    for nq in nqs {
                        set_number_of_fermions(nq);
                        for mq in mqs {
                            let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);
                            for mu in mus {
                                for om in oms {
                                    for p in ps {
                                        assert_equal(
                                            propagator_l_zero_temp_landau(om, p, mg, mq, mu, f0),
                                            propagator_l_zero_temp_landau_w_field_config(
                                                om, p, mu, f0, &config,
                                            ),
                                        )
                                    }
                                }
                            }
                        }
                    }
                }
                set_number_of_colors(3);
                set_number_of_fermions(1);
            }

            #[test]
            #[serial]
            fn test_prop_t_zero_temp_w_field_config() {
                use qcd_sme::qcd::thermal::gluon::{
                    propagator_t_zero_temp_landau, propagator_t_zero_temp_landau_w_field_config,
                };
                use qcd_sme::qcd::FieldConfig;
                use qcd_sme::ym::gluon::propagator_landau as ym_propagator_t_zero_temp_landau;

                let mg = 0.656;
                let f0 = -0.876;

                let ncs = [2, 3];
                let nqs = [1, 2];
                let mqs = [0.1, 0.5];
                let mus = [0., 0.7];
                let oms = [0.37, 1.13];
                let ps = [0.41, 0.94];

                let mg2 = mg * mg;

                for nc in ncs {
                    set_number_of_colors(nc);
                    let config = FieldConfig::new(nc, mg, vec![]);
                    for om in oms {
                        for p in ps {
                            let s = (om * om + p * p) / mg2;
                            assert_equal(
                                ym_propagator_t_zero_temp_landau(s, f0) / mg2,
                                propagator_t_zero_temp_landau_w_field_config(
                                    om, p, 0., f0, &config,
                                )
                                .re,
                            )
                        }
                    }
                }

                for nc in ncs {
                    set_number_of_colors(nc);
                    for nq in nqs {
                        set_number_of_fermions(nq);
                        for mq in mqs {
                            let config = FieldConfig::new(nc, mg, vec![(nq, mq)]);
                            for mu in mus {
                                for om in oms {
                                    for p in ps {
                                        assert_equal(
                                            propagator_t_zero_temp_landau(om, p, mg, mq, mu, f0),
                                            propagator_t_zero_temp_landau_w_field_config(
                                                om, p, mu, f0, &config,
                                            ),
                                        )
                                    }
                                }
                            }
                        }
                    }
                }
                set_number_of_colors(3);
                set_number_of_fermions(1);
            }
        }
    }
}
