# Calculate Mean absolute deviation (error) and Root mean square deviation (error)
# returning only the metric of choice (defaulting to MAD).

resid.metric <- function(ffr.fit, metric = c('MAD', 'RMSD')){
  
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
    
    # absolute deviation
    AD <- abs(Nconsumed - ffr.fit@data$killed)
      
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

    # absoluate deviations
    AD <- c(
      abs(Nconsumed[,1] - ffr.fit@data$Ni_consumed),
      abs(Nconsumed[,2] - ffr.fit@data$Nj_consumed)
    )
    
    # squared deviations
    SD <- c(
      (Nconsumed[,1] - ffr.fit@data$Ni_consumed)^2,
      (Nconsumed[,2] - ffr.fit@data$Nj_consumed)^2
    )
  }
  
  # Mean absolute deviation
  MAD <- mean(AD)
  
  # Root mean square deviation
  RMSD <- sqrt(mean(SD))
  
  # determine which metric to return (defaults to first, i.e. MAE)
  metric <- match.arg(metric)
  
  if(metric == 'MAD'){
    return(MAD)
  }
  
  if(metric == 'RMSD'){
    return(RMSD)
  }
}