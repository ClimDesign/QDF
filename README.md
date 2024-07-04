# Flexible and consistent flood-duration-frequency (QDF) models

Authors: [Danielle M. Barna](https://scholar.google.com/citations?hl=no&user=homV8wQAAAAJ), Kolbjørn Engeland, Thordis Thorarinsdottir, Chong-Yu Xu

Publications: [Flexible and consistent Flood–Duration–Frequency modeling: A Bayesian approach](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)

## Overview
Design values estimate the relationship between a flood’s return level (magnitude) and return period (frequency). Flood-duration-frequency (QDF) models are a class of statistical models that provide design value estimates at multiple durations by simulateneously estimating several durations and quantiles at once under consistency constraints. They are are analogous to intensity-duration-frequency (IDF) models used in precipitation modeling. 

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
Additionally, [Stan](https://mc-stan.org/) (version 2.25 or later) was used to generate initial values for some parameters to improve the efficiency of the MCMC sampler. See (this link) for details. 

### Data structure
QDF analysis requires generating several sets of annual maxima from observed data at a single gauging station. Each set of annual maxima corresponds to a different duration. The model then finds a relationship between these sets of annual maxima. See [Barna et al., 2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34) for a discussion of how to process data for QDF analysis.  

The MCMC sampler expects data in the form of a 2 x (n*d) matrix, where n = number of years of observed data at a station and d = the number of durations to be simultaneously modeled. The first row contains the data values and the second row indicates the duration that the data value corresponds to. The data should be ordered from smallest duration to largest. 

(example? change to Rmd)

### Running the MCMC sampler
For the original and extended QDF models Bayesian inference is performed using a Metropolis-Within-Gibbs algorithm. The mixture model relies on a reversible jump MCMC algorithm, similar to the reversible jump methodology for normal mixtures described in [Richardson and Green (1997)](https://academic.oup.com/jrsssb/article/59/4/731/7083042). See [Barna et al., 2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34), section 3.3 for details.

The three models are called using the commands `javelle()`, `extendedQDF()`, and `reversiblejumpQDF()`. 

Time? output?

## 
