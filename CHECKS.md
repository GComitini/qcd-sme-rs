<!-- markdownlint-disable MD024 -->
# Checks

## Generic

- Low-level functions (logarithmic and rational vacuum functions, thermal integrands, etc.)
  - [x] Numerical tests using independent Python code (module-level `tests` submodules and [`tests/numerics.py`](tests/numerics.py))

## Pure Yang-Mills theory

### Vacuum

- Ghost propagator
- Gluon propagator
  - [x] Values, Landau gauge ([`examples/qcd_gluon_propagator_compare.rs`](examples/qcd_gluon_propagator_compare.rs), [arXiv:1605.07357](https://arxiv.org/pdf/1605.07357.pdf))
  - [x] Poles, Landau gauge ([`examples/ym_gluon_poles.rs`](examples/ym_gluon_poles.rs), [arXiv:1806.08397](https://arxiv.org/abs/1806.08397))

### Finite temperature

- Ghost propagator
- Gluon propagator
  - [x] Temperature-dependent poles ([`examples/ym_gluon_poles.rs`](examples/ym_gluon_poles.rs), [arXiv:2101.08341](https://arxiv.org/abs/2101.08341))
  - [x] Values ([`examples/ym_gluon_propagator.rs`](examples/ym_gluon_propagator.rs), [arXiv:2101.08341](https://arxiv.org/abs/2101.08341))

## Full QCD

### Vacuum

- Ghost propagator
- Gluon propagator
  - [x] Values, Landau gauge ([`examples/qcd_gluon_propagator_compare.rs`](examples/qcd_gluon_propagator_compare.rs), [arXiv:1605.07357](https://arxiv.org/pdf/1605.07357.pdf))

### Finite temperature/chemical potential

- Ghost propagator
- Gluon propagator
  - [x] Relative coefficients of gluon and quark thermal contribution: asymptotic (high temperature/large chemical potential) HTL-like scaling ([`examples/qcd_gluon_propagator_asymptotics.rs`](examples/qcd_gluon_propagator_asymptotics.rs))
  - [x] Zero-momentum limit equalness of polarization projections ([`examples/qcd_gluon_polarization.rs`](examples/qcd_gluon_polarization.rs))