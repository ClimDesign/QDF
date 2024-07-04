# QDF

Authors: [Danielle M. Barna](https://scholar.google.com/citations?hl=no&user=homV8wQAAAAJ), Kolbjørn Engeland, Thordis Thorarinsdottir, Chong-Yu Xu

## Overview
Design values estimate the relationship between a flood’s return level (magnitude) and return period (frequency). Flood-duration-frequency (QDF) models are a class of statistical models that provide design value estimates at multiple durations by simulateneously estimating several durations and quantiles at once under consistency constraints. They are are analogous to intensity-duration-frequency (IDF) models used in precipitation modeling. 

This respository contains routines to fit three different QDF models:
- The original QDF model from [Javelle et al., 2002](https://www.sciencedirect.com/science/article/pii/S0022169401005777)
- An extended "multiscaling" QDF model from [Barna et al., 2023](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34)
- A mixture model that combines elements of both the original and extended models

<img src="https://github.com/ClimDesign/QDF/assets/49793254/ee515bf8-0c44-4436-88e0-44b1ea58c726" width="600">

Return level plots from a synthetic data set showing output from (i) the original "simple scaling" QDF model and (ii) output from a extended "multiscaling" QDF model that relaxes the consistency constraints to allow for more realistic behavior. Figure adapted from Barna et al., (2023).

## Publications
[Flexible and consistent Flood–Duration–Frequency modeling: A Bayesian approach](https://www.sciencedirect.com/science/article/pii/S0022169423003906#b34).

## How to fit the models
The following dependencies are needed to fit the models:
```
library(extRemes)
library(EnvStats)
library(truncnorm)

source("mcmc_sampler.R")
```
Additionally, the efficiency of the mcmc sampler is greatly improved if local 
