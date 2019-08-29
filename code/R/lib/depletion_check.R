depletion.check<-function(ffr.fit, fracEaten=0.9, fracRepl=0.1){
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
  n.depleted <- sum(eaten/initial > fracEaten)
  frac.depleted <- n.depleted/length(initial)
  depleted <- frac.depleted>fracRepl
  return(depleted)
}


# depletion.check(ffr.fits[[1]])
# lapply(ffr.fits, depletion.check)
