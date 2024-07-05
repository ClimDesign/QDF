# Flexible and consistent flood-duration-frequency (QDF) models

Authors: [Danielle M. Barna](https://scholar.google.com/citations?hl=no&user=homV8wQAAAAJ), Kolbjørn Engeland, Thordis Thorarinsdottir, Chong-Yu Xu

Publications: [Flexible and consistent Flood–Duration–Frequency modeling: A Bayesian approach](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)

## Overview
Design values estimate the relationship between a flood’s return level (magnitude) and return period (frequency). Flood-duration-frequency (QDF) models are a class of extreme value models that provide design value estimates at multiple durations by simulateneously estimating several durations and quantiles at once under consistency constraints.  Typically the underlying extreme value distribution is assumed to be the generalized extreme value (GEV) distribution. QDF models are a type of *dependent GEV*, or d-GEV model, and are analogous to intensity–duration–frequency (IDF) models for precipitation.

This respository contains routines to fit three different QDF models:
- The original QDF model from [Javelle et al., 2002](https://www.sciencedirect.com/science/article/pii/S0022169401005777)
- An extended "multiscaling" QDF model from [Barna et al., 2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)
- A mixture model that combines elements of both the original and extended models

<img src="https://github.com/ClimDesign/QDF/assets/49793254/ee515bf8-0c44-4436-88e0-44b1ea58c726" width="600">

Example output: return level plots from a synthetic data set showing output from (i) the original "simple scaling" QDF model and (ii) output from a extended multiscaling QDF model where consistency constraints are relaxed to allow for more realistic behavior. Figure adapted from Barna et al., (2023).

## How to fit the models
The following dependencies are needed to fit the models:
```
library(extRemes)
library(EnvStats)
library(truncnorm)

source("mcmc_sampler.R")
```
The models are fit in a Bayesian framework and estimation relies on Markov Chain Monte Carlo (MCMC) sampling.

[Stan](https://mc-stan.org/) (version 2.25 or later) was used to generate initial values for some parameters to improve the efficiency of the MCMC sampler. Details on using Stan to generate initial values for QDF analysis are available at (this link).

### Data structure
QDF analysis requires generating several sets of annual maxima from observed data at a single gauging station. Each set of annual maxima corresponds to a different duration. The model then finds a relationship between these generated sets of annual maxima. See [Barna et al., (2023)](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34) for a discussion of how to process data for QDF analysis.  

The MCMC sampler expects data in the form of a 2 x (n*d) matrix, where n = number of years of observed data at a station and d = the number of durations to be simultaneously modeled. The first row contains the data values and the second row indicates the duration that the data value corresponds to. The data should be ordered from smallest duration to largest. 

(example? change to Rmd)

### Running the MCMC sampler
For the original and extended QDF models Bayesian inference is performed using a Metropolis-Within-Gibbs algorithm. The mixture model relies on a reversible jump MCMC algorithm, similar to the reversible jump methodology for normal mixtures described in [Richardson and Green (1997)](https://academic.oup.com/jrsssb/article/59/4/731/7083042). The models the algorithm jumps between are the original and extended QDF models. See [Barna et al., 2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34), section 3.3 for details.

The three models are called using the commands `javelle()`, `extendedQDF()`, and `reversiblejumpQDF()`. The inputs are
- `data` - the data in matrix format
- `startpoint` - a named vector containing startpoints for the model parameters. 
- `tuning` - a named vector specifying the width of the proposal distributions for the model parameters.
- `iter` - the number of iterations in the MCMC sampler.
- `ss` - the number of iterations at which we subsample. If included, the sampler stores every ss-th iteration.

The `reversiblejumpQDF()` model additionally takes 
- `innerloop` - the number of times to repeat sampling before proposing to jump between models.

For example, if we wanted to run the `reversiblejumpQDF()` model for 2.5*10^6 iterations, storing every 5th iteration, where we repeat sampling 10 times before proposing a jump between models, we would call:
```
outputRJ <- reversiblejumpQDF(data, startpoint, tuning, iter = 2.5*10^6, ss = 5, innerloop = 10)
```

The routines presented were developed for targeted analysis of twelve gauging stations, allowing individual generation of start points and tuning of proposal distributions. These start points and tuning parameters are stored in the repository. 

The model output is stored as a list, where the first list item stores the log prior and log likelihood values at each iteration and the second list item stores the parameter values at each iteration. The `javelle()` and `extendedQDF()` models also store the acceptance rates for each parameter as a third list item.

## Working with the model output
The models rely on a quantile-based reparameterization of the GEV where the location parameter is replaced with the median. See, e.g., [Castro-Camilo et al., (2022)](https://onlinelibrary.wiley.com/doi/full/10.1002/env.2742). 

To convert between this reparameterization and the more common [location-scale parametrization](https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution) of the GEV, use `this`. (rowwise conversion? mixture model?

To plot return level plots, use `this`. 

## References
