# Deprecated function.  Use resid.metric() in resid.metrics.R instead.

# Root mean square deviation
RMSD <- function(ffr.fit){

  # grab model parameters
  params <- coef(ffr.fit)

  # grab model type
  modeltype <- ffr.fit@data$modeltype

  # check if a 1pred1prey fit
  if("initial" %in% names(ffr.fit@data)){
    # perform parameter transformation
    set_params(params, modeltype)

    if(modeltype %in% c(
        'Holling.I',
        'Holling.II',
        'Beddington.DeAngelis',
        'Crowley.Martin',
        'Stouffer.Novak.I',
        'Stouffer.Novak.II',
        'Stouffer.Novak.III'
        )
    ){
      # predicted consumption
      Nconsumed <- holling.like.1pred.1prey.predict(
        params,
        modeltype=ffr.fit@data$modeltype,
        initial=ffr.fit@data$initial,
        predators=ffr.fit@data$predators,
        replacement=ffr.fit@data$replacement,
        Pminus1=ffr.fit@data$Pminus1,
        time=ffr.fit@data$time
      )
    }else{
      # predicted consumption
      Nconsumed <- ratio.like.1pred.1prey.predict(
        params,
        modeltype=ffr.fit@data$modeltype,
        initial=ffr.fit@data$initial,
        predators=ffr.fit@data$predators,
        replacement=ffr.fit@data$replacement,
        time=ffr.fit@data$time
      )
    }

    # squared deviation
    SD <- (Nconsumed - ffr.fit@data$killed)^2
  }else{
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

    # squared deviations
    SD <- c(
      (Nconsumed[,1] - ffr.fit@data$Ni_consumed)^2,
      (Nconsumed[,2] - ffr.fit@data$Nj_consumed)^2
    )
  }

  # Root mean square deviation
  RMSD <- sqrt(mean(SD))

  return(RMSD)
}