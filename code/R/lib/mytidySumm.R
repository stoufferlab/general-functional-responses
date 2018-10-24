# DEBUG: we should probably move these functions elsewhere so this file only calls functions/runs analysis
# library(devtools)
# devtools::install_github("bbolker/broom")
library(broom) # for tidy()

mytidy <- function(fit){ 
  if(typeof(fit)=='S4'){
    tfit <- tidy(fit)
    terms <- tfit$term
    tfit$term <- NULL
    out <- matrix(as.numeric(unlist(tfit)), nrow=nrow(tfit), ncol=ncol(tfit), byrow=FALSE)
    rownames(out) <- terms
    colnames(out) <- colnames(tfit)
    return(out)
  }else{
    return(NA)}
}

make.array <- function(ffr.fit, boot.reps){
  # if(typeof(ffr.fit)=='logical'){
  #   return(ffr.fit)
  # }
  if(typeof(ffr.fit)=='S4'){
    ffr.fit <- mytidy(ffr.fit)
  }
  out <- array(NA, dim=c(dim(ffr.fit), boot.reps))
  dimnames(out) <- dimnames(ffr.fit)
  return(out)
}

summarize.boots <- function(x){
  c(mean=mean(x, na.rm=TRUE), quantile(x,c(0.025,0.16,0.5,0.84,0.975),na.rm=TRUE), n=sum(!is.na(x)))
}
