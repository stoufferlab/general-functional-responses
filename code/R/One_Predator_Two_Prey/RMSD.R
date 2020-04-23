# Root mean square deviation (error)
# Error technically refers to out-of-sample
# Deviation technically refers to within-sample, hence use RMSD in manuscript
RMSD <- function(ffr.fit){

  # scrape out the parameters
  params <- coef(ffr.fit)

  # generate the predicted values
  Nconsumed <- holling.like.1pred.2prey.predict(
    params,
    Ni=ffr.fit@data$Ni,
    Nj=ffr.fit@data$Nj,
    Npredators=ffr.fit@data$Npredators,
    replacement=ffr.fit@data$replacement,
    modeltype=ffr.fit@data$modeltype,
    phi.transform=ifelse(
      "phi.transform" %in% names(ffr.fit@data),
      ffr.fit@data$phi.transform,
      NULL
    ),
    time=ffr.fit@data$time
  )

  # Root mean square deviation
  RMSD <- sqrt(
    mean(c(
      (Nconsumed[,1]-ffr.fit@data$Ni_consumed)^2,
      (Nconsumed[,2]-ffr.fit@data$Nj_consumed)^2
    ))
  )

  return(RMSD)
}