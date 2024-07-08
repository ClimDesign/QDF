A quantile-based implementation of the GEV distribution in Stan
================

## Overview

We present a user-defined implementation of the generalized extreme
value (GEV) distribution in the probabilistic programming language
[Stan](https://mc-stan.org/).

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
GEV. See, e.g., [Castro-Camilo et al.,
(2022)](https://onlinelibrary.wiley.com/doi/full/10.1002/env.2742).

### Prerequisites

Requires [Stan version
2.25](https://mc-stan.org/rstan/reference/stan_version.html) or greater.

## Implementation details

In its standard location-scale parameterization the GEV has CDF

$$
G(y) = \text{exp}(-[1 + \xi(\frac{y-\mu}{\sigma})]^{-1/\xi}).
$$

The support of the distribution depends on all three parameters and the
data:

$$
\tag{1}
\{ y : 1 + \xi(y-\mu)/\sigma > 0 \}
$$

where $-\infty < \mu < \infty$, $\sigma > 0$, and
$-\infty < \xi < \infty$.

When fitting the GEV with a boilerplate MCMC (for example, a
metropolis-hasting algorithm) we would enforce this condition by simply
setting the log likelihood equal to $-\infty$ when the sampler proposes
a value outside the support.

However, Stan needs parameter bounds defined at initialization, so we
can’t just enforce any support conditions in the likelihood.
Furthermore, the parameters need to be continuously differentiable at
the bounds, so we can’t have bounds in the form of conditional
statements on other parameter values (no if-then statements).

This is typically handled (see Aki Vehtari’s
[work](https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html)
with the Generalized Pareto distribution) by forcing the burden of
support onto a single parameter such that the bounds for this one
parameter are a continuously differentiable function of the other
parameters (and, in the case of the GEV, the maximum and minimum data
points as well). In the GEV distribution it makes sense to have the
$\xi$ parameter be the one that accommodates the support.

### How the sign of $\xi$ controls the support of the GEV distribution

If $\xi$ has bounds dependent on the other two parameters we can think
of the sign of $\xi$ as being determined by the sign of the quantity
$$(y-\mu)/\sigma.$$ Assuming all $y$ are positive, we can consider a toy
scenario where $\mu$ is guaranteed to fall between min($y$) and
max($y$):

|                                                                                                                            |
|----------------------------------------------------------------------------------------------------------------------------|
| **If $(y-\mu)/\sigma > 0$ then $\xi < 0$**                                                                                 |
| Then the $y$ we care about is $y_{max}$ and $\xi$ needs a lower bound such that $$\xi > \frac{-1}{(y_{max}-\mu)/\sigma}.$$ |

|                                                                                                                             |
|-----------------------------------------------------------------------------------------------------------------------------|
| **If $(y-\mu)/\sigma < 0$ then $\xi > 0$**                                                                                  |
| Then the $y$ we care about is $y_{min}$ and $\xi$ needs an upper bound such that $$\xi < \frac{-1}{(y_{min}-\mu)/\sigma}.$$ |

This works because we can trap $\xi$ between two bounds, one always
dependent on $y_{max}$ and the other always dependent on $y_{min}$. But
this logic fails when applied to the GEV distribution in its standard
location-scale parameterization since $\mu$ can fall outside the range
of the data. This means we cannot definitively tie $y_{max}$ to the
lower bound and $y_{min}$ to the upper bound, respectively, since we
would need to know the value of $\mu$ ahead of time to know which order
statistic we care about in the context of Eqn (1).

### Enforcing support via a quantile-based reparameterization of the GEV

The trick is to reparameterize the GEV such that we can make claims
about one of the parameters that allow us to reliably tie the order
statistics ($y_{min}$,$y_{max}$) to the bounds on $\xi$.

The relationship between the location parameter, $\mu$, and the median,
$\eta$, is given as $$
\begin{equation}
    \eta = 
    \begin{cases}
    \mu + \sigma \frac{\text{log}(2)^{-\xi}-1}{\xi} & \text{if}\ \xi \neq 0 \\
    \mu - \text{log}\left(\text{log}(2)\right) & \text{if}\ \xi = 0.
    \end{cases}
\end{equation}
$$
