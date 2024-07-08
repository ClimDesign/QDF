A quantile-based implementation of the GEV distribution in Stan
================

## Overview

The functions in `ifgev.stan` provide a user-defined implementation of
the generalized extreme value (GEV) distribution in the probabilistic
programming language [Stan](https://mc-stan.org/).

Stan provides full Bayesian statistical inference using a powerful
Hamiltonian Monte Carlo (HMC) algorithm, provided the parameter bounds
in the model are properly defined [(Gelman et al.,
2015)](https://journals.sagepub.com/doi/abs/10.3102/1076998615606113).
Stan allows for [very fast sampling](https://arxiv.org/abs/1206.1901),
but, most importantly, provides a [range of useful
diagnostics:](https://mc-stan.org/docs/stan-users-guide/) using HMC to
explore the target distribution means failures in geometric ergodicity
manifest in distinct behaviors that can be developed into diagnostic
tools. Accurate diagnostics are particularly appealing in an extreme
value setting given the natural sparsity of extreme value data.

The implementation relies on a quantile-based reparameterization of the
GEV. Parameter bounds are worked out in [this Rpubs
document](https://rpubs.com/dbarna/brmsgev) which is associated with
[this `brms` issue
report.](https://github.com/paul-buerkner/brms/issues/1345)

## Prerequisites

Requires [Stan version
2.25](https://mc-stan.org/rstan/reference/stan_version.html) or greater.
