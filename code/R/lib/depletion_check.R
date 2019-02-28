depletion.check<-function(ffr.fit, cutoff=0.8){
  # only bother if non-replacement
  if(!ffr.fit$study.info$replacement){
    return(0)
  }else
    if(!is.null(ffr.fit$study.info$data.Nconsumed)){
      eaten <- ffr.fit$study.info$data.Nconsumed
    }else{
      eaten <- ffr.fit$study.info$data.Nconsumed.mean
    }
  initial <- ffr.fit$study.info$data.Nprey
  n.depleted <- sum(eaten/initial > cutoff)
  frac.depleted <- n.depleted/length(initial)
  return(frac.depleted)
}


# depletion.check(ffr.fits[[1]])
# lapply(ffr.fits, depletion.check)
