# Root mean square error
RMSE <- function(d,
                 ffr.fit,
                 study.info,
                 model=c('Holling.I',
                         'Holling.II',
                         'Beddington.DeAngelis',
                         'Crowley.Martin',
                         'Stouffer.Novak.I',
                         'Stouffer.Novak.II',
                         'Stouffer.Novak.III',
                         'Ratio',
                         'Arditi.Ginzburg',
                         'Hassell.Varley',
                         'Arditi.Akcakaya')){

  params <- coef(ffr.fit)
  set_params(params)
  model <- match.arg(model)
  
  # expected number consumed given data and parameters
  if(model %in% c('Holling.I',
                  'Holling.II',
                  'Beddington.DeAngelis',
                  'Crowley.Martin',
                  'Stouffer.Novak.I',
                  'Stouffer.Novak.II',
                  'Stouffer.Novak.III')){
    Nconsumed.predicted <- holling.like.1pred.1prey(N0=d$Nprey, 
                                                    a=attack, 
                                                    h=handling, 
                                                    c=interference, 
                                                    phi_numer=phi_numer, 
                                                    phi_denom=phi_denom, 
                                                    P=d$Npredator, 
                                                    T=d$Time, 
                                                    replacement=study.info$replacement, 
                                                    Pminus1=study.info$Pminus1)
  }else if(model %in% c('Ratio',
                        'Arditi.Ginzburg',
                        'Hassell.Varley',
                        'Arditi.Akcakaya')){
    Nconsumed.predicted <- ratio.like.1pred.1prey(N0=d$Nprey,
                                                  a=attack, 
                                                  h=handling,
                                                  m=exponent, 
                                                  P=d$Npredator, 
                                                  T=d$Time, 
                                                  replacement=study.info$replacement)
  } else {stop("Incorrect model specification in RMSE()")}
  
  # Clean up
  rm(attack, handling, interference, phi_numer, phi_denom, exponent, envir = .GlobalEnv)

  # Root mean square error
  RMSE <- sqrt(mean((d$Nconsumed-Nconsumed.predicted)^2))

  return(RMSE)
}