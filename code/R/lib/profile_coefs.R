
# useful for an implementation to get SEs from non-positive-definite Hessians
library(HelpersMG)

# Function to performing interval profiling on estimated model coefficients for subsequent plotting. Intervals are determined on estimated scale and must be backtransformed to plot on natural scale.  Script applying function it to all datasets follows.
profile_coefs <- function(
  ffr.fits,
  model=c('Holling.I',
          'Holling.II',
          'Beddington.DeAngelis',
          'Crowley.Martin',
          'Stouffer.Novak.I',
          'Ratio',
          'Arditi.Ginzburg',
          'Hassell.Varley',
          'Arditi.Akcakaya',
          'Arditi.Akcakaya.Method.2'),
  point.est=c('median','mean'),
  printWarnings = TRUE,
  which.pars=NULL,
  ...
){
  
  model <- match.arg(model)
  point.est <- match.arg(point.est)
  point.est <- ifelse(point.est=='median', '50%', point.est)
 
  require(progress)
  pb <- progress_bar$new(
    format = "Profiling intervals for each data set [:bar] :percent eta: :eta",
    total = length(ffr.fits),
    show_after = 0,
    force = TRUE,
    clear = FALSE
  )
  
  if(printWarnings){options(warn=1)}else{options(warn=0)}
    
    profiles <- list()
    for(i in 1:length(ffr.fits)){

      x <- ffr.fits[[i]]
      print(x$study.info$datasetName)
      
      if(is.null(which.pars)){
        pars <- colnames(x$estimates[[model]])
      }else{
        pars <- which.pars
      }

      # extract point estimate (the median estimate is easy to determine regardless of the type of data so should be default)
      est <- x$estimates[[model]][point.est, pars, "estimate"]
      
      # cheeky lower and upper bounds to start
      lb <- est
      ub <- est

      # make sure we have a well-behaved Hessian to use to jumpstart profiling
      helper <- try(HelpersMG::SEfromHessian(x$fits[[model]]@details$hessian, hessian=TRUE))
      if(!inherits(helper, "try-error")){
        std.err <- helper$SE
        x$fits[[model]]@details$hessian <- helper$hessian
        x$fits[[model]]@vcov <- solve(helper$hessian)
      }else{
        std.err <- summary(x$fits[[model]])@coef[pars,"Std. Error"]
        std.err[is.na(std.err)] <- 0.01
      }

      # if we did not bootstrap then try (1) profile or (2) approximate	
      if(x$estimates[[model]]["n",1,1][1] == 1){
        
        # (1) estimate the profile confidence interval
        if(model!='Arditi.Akcakaya.Method.2'){
          # message('attempt 1')
          cf <- try(
            confint(
              profile(
                x$fits[[model]],
                which=pars,
                tol.newmin=Inf,
                try_harder=TRUE,
                std.err=std.err,
                ...
              ),
              level=0.68
            )
          )
          # # sometimes profiling fails because it needs specification of the step sizes to use
          # if(inherits(cf, "try-error") || any(is.na(cf))){
          #   # message('attempt 2')
          #   std.err <- summary(x$fits[[model]])@coef[pars,"Std. Error"]
          #   std.err[is.na(std.err)] <- 0.01
          #   cf <- try(
          #     confint(
          #       profile(
          #         x$fits[[model]],
          #         which=pars,
          #         tol.newmin=Inf,
          #         try_harder=TRUE,
          #         std.err=std.err,
          #         # maxsteps=1E6,
          #         ...
          #       ),
          #       level=0.68
          #     )
          #   )
          # }
          # and one last hail mary
          if(inherits(cf, "try-error") || any(is.na(cf))){
            # I cannot for the life of me remember where I got this idea from, but sometimes it works...
            x$fits[[model]]@call$start <- x$fits[[model]]@call$data$sbplx.start
            x$fits[[model]]@call$control <- NULL

            cf <- try(
              confint(
                profile(
                  x$fits[[model]],
                  which=pars,
                  tol.newmin=Inf,
                  try_harder=TRUE,
                  std.err=std.err,
                  ...
                ),
                level=0.68
              )
            )
          }
          # if any are still NA treat as if we got an error
          if(!all(is.finite(cf))){
            cf <- TRUE
            class(cf) <- 'try-error'
          }
        }else{ # if AA method 2 model
          cf <- TRUE
          class(cf) <- 'try-error'
        }
        
        # if profiling code was successful
        if(!inherits(cf, "try-error")){
          method <- 'profile'

          if(length(pars)==1){
            lb <- cf[1]
            ub <- cf[2]
          }else{
            lb <- cf[,1]
            ub <- cf[,2]
          }
          
        }else{
          # (2) if profiling is unsuccessful or AA2 method was used then assume quadratic approximation
          method <- 'quadratic'
          
          # get the SEs directly from the model output
          if(model=='Arditi.Akcakaya.Method.2'){
            se <- x$estimates$Arditi.Akcakaya.Method.2[,,'std.error'][point.est,]
          }else{
            se <- coef(summary(x$fits[[model]]))[pars,"Std. Error"]
          }

          lb <- est - se
          ub <- est + se
        }
      }else{# (3) if we bootstrapped then use the quantiles
        method <- 'bootstrap'

        # use central interval equivalent to one SD as bounds
        lb <- x$estimates[[model]]["16%",pars,"estimate"]
        ub <- x$estimates[[model]]["84%",pars,"estimate"]
      }

      cf <- as.data.frame(cbind(lb,est,ub))
      rownames(cf) <- pars

      profiles[[i]] <- list(
        study.info = x$study.info,
        model = model,
        estimates = x$estimates[model],
        profile=list(
          cf = cf,
          method = method,
          point.est = point.est
        )
      )
      pb$tick()
    }

  options(warn=0)
  return(profiles)
}
