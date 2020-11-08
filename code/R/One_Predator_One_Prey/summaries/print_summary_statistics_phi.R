
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

# Read in the AIC summary table
stat <- '50%' # use "50%" or "mean"

# generate a table of AIC model ranks
source('make_AIC_table_phi.R')

# statistics about when the phi model ranks first or ties for first based on AIC
message("AIC model ranks")
message(paste(
  "phi itself:",
  sum(rnkAICs[,"G"]==1),
  round(sum(rnkAICs[,"G"]==1)/nrow(rnkAICs)*100),
  "%"
))
message(paste(
  "phi tied:",
  sum(dAICs[,"G"]<=2)-sum(rnkAICs[,"G"]==1),
  round((sum(dAICs[,"G"]<=2)-sum(rnkAICs[,"G"]==1))/nrow(rnkAICs)*100),
  "%"
))
message()

# generate a table of AICc model ranks
source('make_AICc_table_phi.R')

# statistics about when the phi model ranks first or ties for first based on AICc
message("AICc model ranks")
message(paste(
  "phi itself:",
  sum(rnkAICcs[,"G"]==1),
  round(sum(rnkAICcs[,"G"]==1)/nrow(rnkAICcs)*100),
  "%"
))
message(paste(
  "phi tied:",
  sum(dAICcs[,"G"]<=2)-sum(rnkAICcs[,"G"]==1),
  round((sum(dAICcs[,"G"]<=2)-sum(rnkAICcs[,"G"]==1))/nrow(rnkAICcs)*100),
  "%"
))
message()

# generate a table of BIC model ranks
source('make_BIC_table_phi.R')

# statistics about when the phi model ranks first or ties for first based on BIC
message("BIC model ranks")
message(paste(
  "phi itself:",
  sum(rnkBICs[,"G"]==1),
  round(sum(rnkBICs[,"G"]==1)/nrow(rnkBICs)*100),
  "%"
))
message(paste(
  "phi tied:",
  sum(dBICs[,"G"]<=2)-sum(rnkBICs[,"G"]==1),
  round((sum(dBICs[,"G"]<=2)-sum(rnkBICs[,"G"]==1))/nrow(rnkBICs)*100),
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