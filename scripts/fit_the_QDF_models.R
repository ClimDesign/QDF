##
##
##
##
##
## Script to run the QDF models and store the output
## -----------------------------------------------------------------------------

library(extRemes)
library(EnvStats)
library(truncnorm)

source("mcmc_sampler.R")

station = "dyrdalsvatn" # station name

load(paste0(station,"_data.rda"))

iter = 2.5*10^6
ss = 10
innerloop = 10

# choose which durations (hours) to fit the model on. See Barna et al., (2023)
# for a discussion of which durations to choose.
durations <- which(stationData[,2] %in% c(1,24,48,72))


# Fit the models ----------------------------------------------------------

# original (simple scaling) QDF model
outputJ <- javelle(stationData[durations,],startpoint,tuning,iter,ss)
saveRDS(outputJ,file=paste0(station,"_1_24_48_72_outputJ.rds"))

# extended (multiscaling) QDF model
outputEx <- extendedQDF(stationData[durations,],startpoint,tuning,iter,ss)
saveRDS(outputEx,file=paste0(station,"_1_24_48_72_outputEx.rds"))

# mixture model. Reversible jump between original and extended models
outputRJ <- reversiblejumpQDF(stationData[durations,],
                              startpoint,tuning,iter,ss,innerloop)
saveRDS(outputRJ,file=paste0(station,"_1_24_48_72_outputRJ.rds"))


