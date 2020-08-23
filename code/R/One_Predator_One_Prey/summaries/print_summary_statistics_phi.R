
# we may not need this but just in case...
# library(bbmle)

# for the bundle_fits function
source('../../lib/plot_coefs.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# which were bootstrapped?
bootstrap <- unlist(lapply(
  ffr.fits,
  function(x){
    x$study.info$bootstrap
  }
))

# print out the summary
message("bootstrapped")
message(paste0("raw: ",sum(!bootstrap)))
message(paste0("not: ",sum(bootstrap)))
message()

# what were the sample sizes?
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
message()

# print out the sample size summary
message("sample sizes")
print(summary(sample.sizes))
message()

# Grab summary of AIC estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"
AIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.I'][[1]][stat]}))
AIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.II'][[1]][stat]}))
AIC.BD <- unlist(lapply(ffr.fits, function(x){ x$AIC['Beddington.DeAngelis'][[1]][stat]}))
AIC.CM <- unlist(lapply(ffr.fits, function(x){ x$AIC['Crowley.Martin'][[1]][stat]}))
AIC.SN1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Stouffer.Novak.I'][[1]][stat]}))

AICs <- data.frame(AIC.H1, AIC.H2, AIC.BD, AIC.CM, AIC.SN1)
colnames(AICs) <- sub('AIC.', '', colnames(AICs))
colnames(AICs)[5] <- "Gen"
# rownames(AICs) <- labels

# make things relative
minAICs <- apply(AICs, 1, min)
dAICs <- AICs - minAICs
# dAICs[dAICs<2] <- 0

# get model ranks
rnkAICs <- t(apply(dAICs, 1, rank, ties.method='first'))
colnames(rnkAICs) <- colnames(AICs)

# statistics about when the phi model ranks first or ties for first
message("model ranks")
message(paste(
  "phi itself:",
  sum(rnkAICs[,"Gen"]==1),
  round(sum(rnkAICs[,"Gen"]==1)/nrow(rnkAICs)*100),
  "%"
))
message(paste(
  "phi tied:",
  sum(dAICs[,"Gen"]<=2)-sum(rnkAICs[,"Gen"]==1),
  round((sum(dAICs[,"Gen"]<=2)-sum(rnkAICs[,"Gen"]==1))/nrow(rnkAICs)*100),
  "%"
))
message()

# load profiled fits to find numbers of datasets with overlap over predefined regions
load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')

message("phi regions")

# count points estimates less than or equal to one
message(paste(
  "point estimate <= 1:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      return(x$profile$cf$est <= 1)
    }
  ))),
  "out of",
  length(ffr.fits)
))

# CI only overlaps beddington deangelis (phi=1)
message(paste(
  "BD only:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(0 < lb && lb <= 1 && 1 <= ub)
    }
  )),
  na.rm=TRUE)
))

# CI only overlaps crowley martin (phi=0)
message(paste(
  "CM only:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(lb < 0 && 0 < ub && ub < 1)
    }
  )),
  na.rm=TRUE)
))

# CI overlaps both crowley martin and beddington deangelis (phi=0 or phi=1)
message(paste(
  "BD and CM:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(lb < 0 && 1 < ub)
    }
  )),
  na.rm=TRUE)
))

# fully below zero
message(paste(
  "phi totally below 0:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(lb < 0 && ub < 0)
    }
  )),
  na.rm=TRUE)
))

# fully between CM and BD models
message(paste(
  "between BD and CM:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(0 < lb && ub < 1)
    }
  )),
  na.rm=TRUE)
))

# fully above one
message(paste(
  "phi totally above 1:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(1 < lb && 1 < ub)
    }
  )),
  na.rm=TRUE)
))

# includes 0
message(paste(
  "phi includes 0:",
  sum(unlist(lapply(
    ffr.fits,
    function(x){
      lb <- ifelse(
        is.finite(x$profile$cf$lb),
        x$profile$cf$lb,
        -Inf
      )
      ub <- ifelse(
        is.finite(x$profile$cf$ub),
        x$profile$cf$ub,
        Inf
      )
      return(lb < 0 && 0 < ub)
    }
  )),
  na.rm=TRUE)
))