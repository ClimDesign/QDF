Flexible and consistent flood-duration-frequency (QDF) models
================

Authors: [Danielle M.
Barna](https://scholar.google.com/citations?hl=no&user=homV8wQAAAAJ),
Kolbjørn Engeland, Thordis Thorarinsdottir, Chong-Yu Xu

Publications: [Flexible and consistent Flood–Duration–Frequency
modeling: A Bayesian
approach](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)

## Overview

Design values estimate the relationship between a flood’s return level
(magnitude) and return period (frequency). Flood-duration-frequency
(QDF) models are a class of extreme value models that provide design
value estimates at multiple durations by simulateneously estimating
several durations and quantiles at once under consistency constraints.
Typically the underlying extreme value distribution is assumed to be the
generalized extreme value (GEV) distribution. QDF models are a type of
*dependent GEV*, or d-GEV model, and are analogous to
intensity–duration–frequency (IDF) models for precipitation.

This respository contains routines to fit three different QDF models:

- The original QDF model from [Javelle et al.,
  2002](https://www.sciencedirect.com/science/article/pii/S0022169401005777)

- An extended “multiscaling” QDF model from [Barna et al.,
  2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)

- A mixture model that combines elements of both the original and
  extended models

<img src="https://github.com/ClimDesign/QDF/assets/49793254/ee515bf8-0c44-4436-88e0-44b1ea58c726" width="600">

Example output: return level plots from a synthetic data set showing
output from (i) the original “simple scaling” QDF model and (ii) output
from a extended multiscaling QDF model where consistency constraints are
relaxed to allow for more realistic behavior. Figure adapted from Barna
et al., (2023).

## How to fit the models

The following dependencies are needed to fit the models:

    library(extRemes)
    library(EnvStats)
    library(truncnorm)

    source("mcmc_sampler.R")

The models are fit in a Bayesian framework and estimation relies on
Markov Chain Monte Carlo (MCMC) sampling.

A custom GEV model defined in [Stan](https://mc-stan.org/) (version 2.25 or later) was used to 
generate initial values for some parameters to improve the efficiency of
the MCMC sampler. The Stan model can be found [here](/Stan).

### Data structure

QDF analysis requires generating several sets of annual maxima from observed data at a single gauging station. Each set of annual maxima corresponds to a different duration. Typically, the duration is defined as the time window over which a specific total flow volume is observed. This focus on total flow volume is useful in engineering design for retention-specific applications, which often require flood volumes for predetermined durations, sometime averaged over different flood events, rather than the multivariate variability of specific flood events. See [Barna et al.,
(2023)](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)
for details.

The MCMC sampler expects data in the form of a (n\*d) x 2 matrix, where
n = number of years of observed data at a station and d = the number of
durations to be simultaneously modeled. The first column contains the
data values and the second column indicates the duration that the data
value corresponds to. The data should be ordered from smallest duration
to largest.

``` r
load("dyrdalsvatn_data.rda")
print(annmax, nrows = 5, digits = 3)
```

    ##     ann max [m^3/s] duration [hours]
    ##               <num>            <num>
    ##  1:            5.77                1
    ##  2:            6.81                1
    ##  3:            7.64                1
    ##  4:            6.68                1
    ##  5:            8.21                1
    ## ---                                 
    ## 92:            5.75               24
    ## 93:            7.92               24
    ## 94:            5.77               24
    ## 95:            6.48               24
    ## 96:            5.42               24

### Running the MCMC sampler

For the original and extended QDF models Bayesian inference is performed
using a Metropolis-Within-Gibbs algorithm. The mixture model relies on a
reversible jump MCMC algorithm, similar to the reversible jump
methodology for normal mixtures described in [Richardson and Green
(1997)](https://academic.oup.com/jrsssb/article/59/4/731/7083042). The
models the algorithm jumps between are the original and extended QDF
models. See [Barna et al.,
2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34),
section 3.3 for details.

The three models are called using the commands `javelle()`,
`extendedQDF()`, and `reversiblejumpQDF()`. The inputs are

- `data` - the matrix of annual maxima and durations

- `startpoint` - a named vector containing startpoints for the model
  parameters.

- `tuning` - a named vector specifying the width of the proposal
  distributions for the model parameters.

- `iter` - the number of iterations in the MCMC sampler.

- `ss` - the number of iterations at which we subsample. If included,
  the sampler stores every ss-th iteration.

The `reversiblejumpQDF()` model additionally takes

- `innerloop` - the number of times to repeat sampling before proposing
  to jump between models.

For example, if we wanted to run the `reversiblejumpQDF()` model for
2.5\*10^6 iterations, storing every 5th iteration, where we repeat
sampling 10 times before proposing a jump between models, we would call:

    outputRJ <- reversiblejumpQDF(annmax, startpoint, tuning, iter = 2.5*10^6, ss = 5, innerloop = 10)

The routines presented were developed for targeted analysis of twelve
gauging stations, allowing individual generation of start points and
tuning of proposal distributions. These start points and tuning
parameters are stored in the repository as part of the `.rda` file for
each station.

The model output is stored as a list, where the first list item stores
the log prior and log likelihood values at each iteration and the second
list item stores the parameter values at each iteration. The `javelle()`
and `extendedQDF()` models also store the acceptance rates for each
parameter as a third list item.

## Working with the model output

The models rely on a quantile-based reparameterization of the GEV where
the location parameter is replaced with the median. See, e.g.,
[Castro-Camilo et al.,
(2022)](https://onlinelibrary.wiley.com/doi/full/10.1002/env.2742).

<DRAFT> To convert between this reparameterization and the more common
[location-scale
parametrization](https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution)
of the GEV, use `this`. (rowwise conversion? mixture model?
