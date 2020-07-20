rm(list = ls())

# we may need a function in this for some reason
source('../../lib/plot_coefs.R')

# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredTwoPrey_fits')

# which datasets were bootstrapped?
bootstrap <- unlist(lapply(
  ffr.fits,
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
  ffr.fits,
  function(x){
    d <- x$study.info
    if("data.Nconsumed1.mean" %in% names(d)){
      return(sum(d$data.n))
    }else{
      return(length(d$data.Nconsumed1))
    }
  }
))
message()

# print out the sample size summary
message("sample sizes")
print(summary(sample.sizes))
message()

# generate the AIC table
source('make_AIC_table.R')

# statistics about when the phi model ranks first or ties for first based on AIC
message("model ranks")
message(paste(
  "phi itself:",
  sum(rnkAICs[,"H2.HH"]==1),
  round(sum(rnkAICs[,"H2.HH"]==1)/nrow(rnkAICs)*100),
  "%"
))
message(paste(
  "phi tied:",
  sum(dAICs[,"H2.HH"]<=2)-sum(rnkAICs[,"H2.HH"]==1),
  round((sum(dAICs[,"H2.HH"]<=2)-sum(rnkAICs[,"H2.HH"]==1))/nrow(rnkAICs)*100),
  "%"
))
message()
