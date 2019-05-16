

# Make set_params() available
sp <- list.files("../../..", "set_params.R", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
source(sp)

# Function to plot observered vs. predicted and calculate pseudo-R2 for each data set and model
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
  Pminus1 <- ffr.fit$study.info$Pminus1
  
  time <- ffr.fit$study.info$data.Time
  # if no times are specified then set to time=1
  if(is.null(time)){
    time <- 1
  }

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
  
  set_params(modeltype, params)
  
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
  legend('bottomright', legend=bquote(R^2==.(round(R2,2))), bty='n',inset=0,cex=0.8)
  title(title,cex=0.5,line=1)
}
