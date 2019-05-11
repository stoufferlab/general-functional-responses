# Function to performing interval profiling on estimated model coefficients for subsequent plotting. Intervals are determined on estimated scale and must be backtransformed to plot on natural scale.  Script applying function it to all datasets follows.
profile_coefs <- function(ffr.fits,
                          model=c('Holling.Type.I',
                                  'Holling.Type.II',
                                  'Beddington.DeAngelis',
                                  'Crowley.Martin',
                                  'Stouffer.Novak.I',
                                  'Ratio',
                                  'Arditi.Ginzburg',
                                  'Hassell.Varley',
                                  'Arditi.Akcakaya',
                                  'Arditi.Akcakaya.Method.2'),
                          point.est=c('median','mean'),
                          ... ){
  
  model <- match.arg(model)
  point.est <- match.arg(point.est)
    point.est <- ifelse(point.est=='median', '50%', point.est)
 
  require(progress)
  pb <- progress_bar$new(
    format = "Profiling intervals for each data set [:bar] :percent eta: :eta",
    total = length(ffr.fits),
    show_after = 0,
    force = TRUE,
    clear = FALSE)
    
  profiles <- list()
  for(i in 1:length(ffr.fits)){
      
    x <- ffr.fits[[i]]
    
    # extract point estimate (the median estimate is easy to determine regardless of the type of data so should be default)
    est <- x$estimates[[model]][point.est, , "estimate"]
    
    # cheeky lower and upper bounds to start
    lb <- est
    ub <- est
    
    # if we did not bootstrap then try (1) profile or (2) approximate	
    if(x$estimates[[model]]["n",1,1][1] == 1){
      
      # (1) estimate the profile confidence interval
      # do so for all model parameters because doing so for focal parameter can cause errors
      if(model!='Arditi.Akcakaya.Method.2'){
        cf <- try(confint(x$fits[[model]], 
                          try_harder=TRUE, 
                          level=0.68, 
                          tol.newmin=Inf, 
                          quietly=TRUE))
      }else{ # if profiling failed in general
        cf <- TRUE
        class(cf) <- 'try-error'
      }
      
      # if profiling code was successful
      if(!inherits(cf, "try-error")){
        method <- 'profile'
       
        lb <- cf[,1]
        ub <- cf[,2]
        
      }else{
        # (2) if profiling is unsuccessful or AA2 method was used then assume quadratic approximation
        method <- 'quadratic'
        
        # get the SEs directly from the model output
        if(model=='Arditi.Akcakaya.Method.2'){
          
          se <- x$estimates$Arditi.Akcakaya.Method.2[,,'std.error'][point.est,]
          
        }else{
          se <- coef(summary(x$fits[[model]]))[,"Std. Error"]
        }
        
        lb <- est - se
        ub <- est + se
        
      }
    }else{# (3) if we bootstrapped then use the quantiles
      method <- 'bootstrap'
  
      # use central interval equivalent to one SD as bounds
      lb <- x$estimates[[model]]["16%",,"estimate"]
      ub <- x$estimates[[model]]["84%",,"estimate"]
      
    }
    
    cf <- cbind(lb,est,ub)
    
    profiles[[i]] <- list(study.info = x$study.info,
                          model = model,
                          estimates = x$estimates[model],
                          profile=list(cf = cf,
                                       method = method,
                                       point.est = point.est))
    pb$tick()
  }
  return(profiles)
}