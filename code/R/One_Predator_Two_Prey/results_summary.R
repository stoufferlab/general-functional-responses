rm(list = ls())

library(bbmle)
# source('../lib/profile_coefs.R')
# for 
source('../lib/plot_coefs.R')
# source('../lib/depletion_check.R') 
# source('../lib/holling_method_one_predator_one_prey.R')
# source('../lib/ratio_method_one_predator_one_prey.R')

# library(RColorBrewer)
# library(Hmisc) # for LaTeX table export
# library(stringr)
# options(xdvicmd='open')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the profiled CIs
load(
  file='../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Stouffer.Novak.I", order.parm="phi_denom")
# ffr.fits <- ffr.fits[fit.order]

# which were bootstrapped?
bootstrap <- unlist(lapply(
  ffr.cfs,
  function(x){
    "data.Nconsumed1.mean" %in% names(x$study.info)
  }
))

# print out the summary
message("bootstrapped")
message(paste0("raw: ",sum(!bootstrap)))
message(paste0("not: ",sum(bootstrap)))

# what were the sample sizes?
sample.sizes <- unlist(lapply(
  ffr.cfs,
  function(x){
    d <- x$study.info
    if("data.Nconsumed1.mean" %in% names(d)){
      return(sum(d$data.n)*2)
    }else{
      return(length(d$data.Nconsumed1)*2)
    }
  }
))

# print out the sample size summary
message("sample sizes")
print(summary(sample.sizes))
