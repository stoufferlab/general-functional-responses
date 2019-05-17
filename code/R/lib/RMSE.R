# Root mean square error
RMSE <- function(d,
                 ffr.fit,
                 study.info,
                 modeltype=c('Holling.type','Ratio.type')){

  params <- coef(ffr.fit)
  set_params(params)
  modeltype <- match.arg(modeltype)
  
  # expected number consumed given data and parameters
  if(modeltype =='Holling.type'){
    Nconsumed.predicted <- holling.like.1pred.1prey(N0=d$Nprey, 
                                                    a=exp(attack), 
                                                    h=exp(handling), 
                                                    c=exp(interference), 
                                                    phi_numer=phi_numer, 
                                                    phi_denom=phi_denom, 
                                                    P=d$Npredator, 
                                                    T=d$Time, 
                                                    replacement=study.info$replacement, 
                                                    Pminus1=study.info$Pminus1)
  }else{
    Nconsumed.predicted <- ratio.like.1pred.1prey(N0=d$Nprey,
                                                  a=exp(attack), 
                                                  h=exp(handling),
                                                  m=exp(exponent), 
                                                  P=d$Npredator, 
                                                  T=d$Time, 
                                                  replacement=study.info$replacement)
  }
  
  # Clean up
  rm(attack, handling, interference, phi_numer, phi_denom, exponent, envir = .GlobalEnv)

  # Root mean square error
  RMSE <- sqrt(mean((d$Nconsumed-Nconsumed.predicted)^2))

  return(RMSE)
}