plot_obsVfit <- function(ffr.fit, 
                         modeltype=c('Holling.Type.I',
                                     'Holling.Type.II',
                                     'Ratio',
                                     'Beddington.DeAngelis',
                                     'Crowley.Martin',
                                     'Arditi.Ginzburg',
                                     'Hassell.Varley',
                                     'Arditi.Akcakaya',
                                     'Stouffer.Novak.I',
                                     'Stouffer.Novak.II',
                                     'Stouffer.Novak.III'),
                         title=c('modeltype','dataset',NULL),
                         ...
){
  
  # Get the data in order
  dataset <- ffr.fit$study.info$datadir
  replacement <- ffr.fit$study.info$replacement
  predators <- ffr.fit$study.info$data.Npredator
  initial <- ffr.fit$study.info$data.Nprey
  time <- ffr.fit$study.info$data.Time
  # if no times are specified then normalize to time=1
  if(is.null(time)){
    time <- 1
  }
  Pminus1 <- ffr.fit$study.info$Pminus1
  
  if(!is.null(ffr.fit$study.info$data.Nconsumed.se)){ # can't use bootstrap to specify since some bootstrapped data don't have SE values and were hence treated as data.
    eaten <- ffr.fit$study.info$data.Nconsumed.mean
    eaten.se <- ffr.fit$study.info$data.Nconsumed.se
  }else{
    eaten <- ffr.fit$study.info$data.Nconsumed
    eaten.se <- NA
  }
  if(title=='modeltype'){title <- modeltype}
  if(title=='dataset'){title <- dataset}
  
  
  # Get the chosen model's predictions
  modeltype <- match.arg(modeltype)
  params <- coef(ffr.fit$fits[[modeltype]])
  
  if(modeltype == "Holling.Type.I"){
    attack <- exp(params[1])
    handling <- 0
    interference <- 0
    phi_numer <- 1
    phi_denom <- 1
  }
  
  if(modeltype == "Holling.Type.II"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    interference <- 0
    phi_numer <- 1
    phi_denom <- 1
  }
  
  if(modeltype == "Ratio"){
    attack <- exp(params[1])
    handling <- 0
    exponent <- 1
  }
  
  if(modeltype == "Arditi.Ginzburg"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    exponent <- 1
  }
  
  if(modeltype == "Hassell.Varley"){
    attack <- exp(params[1])
    handling <- 0
    exponent <- exp(params[2])
  }
  
  if(modeltype == "Arditi.Akcakaya"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    exponent <- exp(params[3])
  }
  
  if(modeltype == "Beddington.DeAngelis"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    interference <- exp(params[3])
    phi_numer <- 1
    phi_denom <- 1
  }
  
  if(modeltype == "Crowley.Martin"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    interference <- exp(params[3])
    phi_numer <- 1
    phi_denom <- 0
  }
  
  if(modeltype == "Stouffer.Novak.I"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    interference <- exp(params[3])
    phi_numer <- 1
    phi_denom <- params[4]
  }
  
  if(modeltype == "Stouffer.Novak.II"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    interference <- exp(params[3])
    phi_numer <- params[4]
    phi_denom <- 1
  }
  
  if(modeltype == "Stouffer.Novak.III"){
    attack <- exp(params[1])
    handling <- exp(params[2])
    interference <- exp(params[3])
    phi_numer <- params[4]
    phi_denom <- params[5]
  }
  
  
  # expected number consumed given data and parameters
  if(modeltype %in% c('Ratio','Arditi.Ginzburg','Hassell.Varley','Arditi.Akcakaya')){
    Nconsumed.predicted <- ratio.like.1pred.1prey(N0=initial,
                                                  a=attack, 
                                                  h=handling,
                                                  m=exponent, 
                                                  P=predators, 
                                                  T=time, 
                                                  replacement=replacement)
  }else{
    Nconsumed.predicted <- holling.like.1pred.1prey(N0=initial, 
                                                    a=attack, 
                                                    h=handling, 
                                                    c=interference, 
                                                    phi_numer=phi_numer, 
                                                    phi_denom=phi_denom, 
                                                    P=predators, 
                                                    T=time, 
                                                    replacement=replacement, 
                                                    Pminus1=Pminus1)
  }
  
  R2 <- cor(eaten, Nconsumed.predicted)^2
  
  # extract LL from (last) fit
  LL <- logLik(ffr.fit$fits[[modeltype]])
  
  rng <- range(c(eaten, Nconsumed.predicted))
  par(pty='s')
  plot(eaten, Nconsumed.predicted,
       ylim=rng,
       xlim=rng,
       ylab='Predicted',
       xlab='Observed',
       type='n')
  if(ffr.fit$study.info$bootstrap){
    arrows(eaten-eaten.se, Nconsumed.predicted, 
           eaten+eaten.se, Nconsumed.predicted, 
           angle=90, length=0.02, code=3)
  }
  abline(0,1)
  points(eaten, Nconsumed.predicted)
  legend('topleft', legend=bquote(R^2==.(round(R2,2))), bty='n',inset=0,cex=0.8)
  # legend('bottomright',legend=bquote(LL==.(round(LL,1))), bty='n',inset=0,cex=0.8)
  title(title,cex=0.5,line=1)
}
