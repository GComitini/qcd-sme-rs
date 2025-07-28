<!-- markdownlint-disable MD024 -->
# Notes

## `qcd_sme`

- The full QCD gluon propagator uses a renormalization which diverges in the $m_{q}\to 0$ limit. To avoid this divergence, the additive renormalization constant $f_{0}$ must be tuned accordingly (e.g. by absorbing in it terms of the form $\ln(m_{q})$).

## `qcd_t_lattice`

- The full QCD gluon propagator fitting functions use the same renormalization as the `qcd_sme` library. When using different conventions (e.g. that in which renormalization is finite in the $m_{q}\to 0$ limit), the additive renormalization constant $f_{0}$ may need to be converted.