A quantile-based implementation of the GEV distribution in Stan
================

In its location-scale parameterization the GEV has CDF

$$
G(y) = \text{exp}\left\{-\left[1 + \xi\left(\frac{y-\mu}{\sigma}\right)\right]^{-1/\xi}\right\}.
$$

$$
G(y) = \text{exp}\{-[1 + \xi(\frac{y-\mu}{\sigma})]^{-1/\xi}\}.
$$

The support of the distribution depends on all three parameters and the
data:

$$
\begin{equation}
\tag{(1)}
\{ y : 1 + \xi(y-\mu)/\sigma > 0 \}
\end{equation}
$$

where $-\infty < \mu < \infty$, $\sigma > 0$, and
$-\infty < \xi < \infty$.

When fitting the GEV with a standard boilerplate MCMC (for example, a
metropolis-hasting algorithm) we would enforce this condition by simply
setting the log likelihood equal to $-\infty$ when the sampler proposes
a value outside the support.

However, once we look at fitting the GEV with more complex methods the
parameter-dependent support starts to become an issue.

In this vignette we look at the gradient-based Hamiltonian Monte Carlo
(HMC) sampler implemented in Stan. Stan needs parameter bounds defined
at initialization, so we can’t just enforce any support conditions in
the likelihood. Furthermore, the parameters need to be continuously
differentiable at the bounds, so we can’t have bounds in the form of
conditional statements on other parameter values (no if-then
statements).

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
$(y-\mu)/\sigma$. Assuming all $y$ are positive, we can consider a toy
scenario where $\mu$ is guaranteed to fall between min($y$) and
max($y$):

<center>

|                                                                                                                            |
|----------------------------------------------------------------------------------------------------------------------------|
| **If $(y-\mu)/\sigma > 0$ then $\xi < 0$**                                                                                 |
| Then the $y$ we care about is $y_{max}$ and $\xi$ needs a lower bound such that $$\xi > \frac{-1}{(y_{max}-\mu)/\sigma}.$$ |

|                                                                                                                             |
|-----------------------------------------------------------------------------------------------------------------------------|
| **If $(y-\mu)/\sigma < 0$ then $\xi > 0$**                                                                                  |
| Then the $y$ we care about is $y_{min}$ and $\xi$ needs an upper bound such that $$\xi < \frac{-1}{(y_{min}-\mu)/\sigma}.$$ |

</center>

This works because we can trap $\xi$ between two bounds, one always
dependent on $y_{max}$ and the other always dependent on $y_{min}$. But
this logic fails when applied to the GEV distribution in its standard
location-scale parameterization since $\mu$ can fall outside the range
of the data. This means we cannot definitively tie $y_{max}$ to the
lower bound and $y_{min}$ to the upper bound, respectively, since we
would need to know the value of $\mu$ ahead of time to know which order
statistic we care about in the context of Eqn (1).
