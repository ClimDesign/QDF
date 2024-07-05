##
##
##
##
##
##
##
## get individual GEV estimates from Stan using annual maximum data
## ----------------------------------------------------------------

library(data.table)
library(rstan)

# set stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# set wd to where stan model lives
setwd("Stan/")

## "annmax" contains the sets of annual maximum from each duration
load("data/dyrdalsvatn_data.rda")

## pick a duration to estimate the Stan model on
duration = 1 # here we pick the 1 hour duration
y <- annmax[which(annmax[,2]==duration),1] 

## ---- fit the model ---------------

# data the Stan model takes is a list where
# N = years of data
# y = observed data points
standata <- list(N = length(y),
                 y = y)

fit <- stan(
  file = 'ifgev.stan',
  data = standata,
  chains = 1,
  warmup = 100,
  iter = 1000,
  cores = 1,
  control = list(adapt_delta = 0.9999),
  seed = 42)

## run a few diagnostics to check the Stan model did okay
monitor(as.array(fit))
check_hmc_diagnostics(fit)

## save the estimated parameters as data.table 'sobj'
sobj <- as.data.table(summary(fit)$summary)
sobj[,param:=c("qind","beta","xi","sigma","mu","lp__")]

## ---- if we want to use Stan output to make     ---------------
## ---- tuning files for fitting the QDF model... ---------------

## startpoints for eta (Qind), beta, xi are the posterior means
## from the Stan fit. Starting values for the Deltas (not available 
## from Stan fit) are chosen as 0.01 and 0.001

startpoint <- matrix(c(sobj[param %in% c("qind","beta","xi"),
                            get("mean")],0.01,0.001),
                     nrow = 1, ncol = 5)
names(startpoint) <- c("Q","B","X","D1","D2")


## the width of the proposal distributions for eta (Qind) is the 
## standard deviations of the posteriors. Standard deviations for the
## Deltas are 0.0001

tuning <- c(0.0001,0.0001,0.0001,sobj[param == "qind",get("sd")])
names(tuning) <- c("Delta1SD","Delta2SD","DeltaJSD","qIndSD")

