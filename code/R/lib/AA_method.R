#######################################################################
# Implementation of 'Method II' of Arditi & Akcakaya 1990.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The method fits a type II functional response to each predator level assuming a common handling time (i.e. it estimates an attack rate for each level), then regresses the log-transformed attack rates on log-transformed predator levels to estimate 'm', weighting each predator level by the inverse of its attack rate estimate's uncertainty (variance).
#######################################################################
require(lamW)
require(bbmle)
require(nloptr)
source('../lib/holling_method_one_predator_one_prey.R') # to fit Holling Type 2
#######################################################################
# ~~~~~~~~~~~~~~~~
# Helper functions
# ~~~~~~~~~~~~~~~~
# Function to determine whether dataset is suitable for Arditi-Akcakaya method.  The method assumes each level of P has variation in N; the data  needs to have at least 2 predator levels as well as at least 3 prey densities for each level (since type II is assumed).
okay4AAmethod<-function(d,minNlevels=3,minPlevels=3){
  tbl <- table(d$Nprey,d$Npredator)
  mns <- apply(tbl>0,2,sum)
  if(all(c(mns>minNlevels, length(mns)>minPlevels)))
    return(TRUE)
  else(return(FALSE))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Negative log likelihood for AAMethod2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AAM.NLL = function(
  params,
  initial,
  killed,
  predators,
  time,
  replacement
){
  
  # P-specific attack rates first and overall handling time last...
  attack <- params[-length(params)]
  handling <- params[length(params)]
  
 # repeat each 'a' for each 'P' level
  MM <- model.matrix(~0+as.factor(predators))
  attack <- MM%*%attack
  
  # we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
  attack <- exp(attack)
  handling <- exp(handling)
  
  # if no times are specified then normalize to time=1
  if(is.null(time)){
    time <- 1
  }
  
  # Expected number consumed given data and parameters
  # The method assumes Holling Type II functional response at each predator abundance level
  # Use holling.like.1pred.1prey() and reduce to Holling type 2
  Nconsumed <- holling.like.1pred.1prey(N0=initial, a=attack, h=handling, c=0, phi_numer=1, phi_denom=1, P=predators, T=time, replacement=replacement, Pminus1=FALSE, overrideTranscendental=TRUE)
  
  # If the parameters are not biologically plausible, neither should be the likelihood
  if(any(Nconsumed < 0 | is.nan(Nconsumed))){
    return(Inf)
  }
  
  # negative log likelihood based on proportion consumed (no replacement)
  if(!replacement){
    nll <- -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE))
  }
  
  # negative log likelihood based on total number consumed (replacement)
  if(replacement){
    nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
  }
  
  return(nll)
}

#######################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main function for Arditi & Akcakaya method
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AAmethod<-function(d,replacement){
    
  nP <- length(unique(d$Npredator))
  xa0 <- log(coef(lm(d$Nconsumed ~ I(d$Nprey*d$Npredator)))[2])
  xh0 <- log(1/quantile(d$Nconsumed,0.9))
  
  fit.AAM.sbplx <- nloptr::sbplx(
    x0=c(rep(xa0,nP),xh0),
    fn=AAM.NLL,
    initial=d$Nprey,
    killed=d$Nconsumed,
    predators=d$Npredator,
    time=d$Time,
    replacement=replacement
  )
  
  if(!is.null(fit.AAM.sbplx$message)){
     fit.AAM.sbplx <- nloptr::sbplx(
      x0=c(rnorm(nP,xa0,0.05),xh0), # add noise starting values
      fn=AAM.NLL,
      initial=d$Nprey,
      killed=d$Nconsumed,
      predators=d$Npredator,
      time=d$Time,
      replacement=replacement,
      control=list(maxeval=2000, xtol_rel=1e-8)
    )
  }
  
  AAM.start<-fit.AAM.sbplx$par
  names(AAM.start) <- parnames(AAM.NLL) <- c(paste0('a',1:nP),'h')
  
  # Refit with mle with nloptr starting estimates 
  # (for convenience to get errors)
  fit.AAM.mle <- bbmle::mle2(
    AAM.NLL,
    start=AAM.start,
    data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, replacement=replacement),
    vecpar=TRUE
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Estimate 'm' as slope of log(a's) ~ log(P's)
  # using "the reciprocal of the variance of y as the weight w."
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ests <- coef(summary(fit.AAM.mle))
  ests.a <- ests[1:nP,1]
  Ps <- unique(attributes(fit.AAM.mle)$data$predators) # grab from fit to ensure Ps are in order corresponding to ests
  
  # Determine weights for each estimate
  w <- 1/diag(vcov(fit.AAM.mle))[1:nP]
  
  # In some cases the variance of the estimates cannot be estimated (resulting either in 0 or NA for the std. error)
  # For these estimates, set weights equal (if none of the SEs can be estimates) or 
  # use the weights of the estimates with the largest variance (lowest weight)
  dg <- diag(vcov(fit.AAM.mle))[1:nP]
  if(all(dg<0 | is.na(dg))){
    dg<-rep(1,length(dg)) # set all weights equal
    w <- 1/dg
    warning('AAMethod2: Unweighted regression used..')
  }
  if(any(dg<0 | is.na(dg))){
    dg[dg<=0|is.na(dg)] <-  max(dg[dg>0|!is.na(dg)],na.rm=T) * 2 # replace twice the largest estimated variance
    w <- 1/dg
    warning('AAMethod2: Twice the largest variance estimate substituted.')
  }

  # est.a estimates are already log-transformed by likelihood function
  fit.AAM.lm <- lm(ests.a ~ log(Ps), weights=w)
  ests.m <- coef(summary(fit.AAM.lm))
  rownames(ests.m) <- c('a0','exponent')
  ests.m[2,1] <- -ests.m[2,1] # negate estimated (negative slope) exponent for interpretation

  # Combine estimates into a single output
  out.ests <- rbind(ests,ests.m)
  colnames(out.ests) <- c('estimate','std.error','statistic','p.value') # for consistency with other nLL fitting output
  
  out <- list(estimates=out.ests, Ps=Ps, fit.a=fit.AAM.mle, fit.m=fit.AAM.lm)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convenience plotting function for Arditi & Akcakaya method
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.AAmethod<-function(AAmethod.out){
  Ps <- AAmethod.out$Ps
  nP <- length(Ps)
  ests.a <- AAmethod.out$estimates[1:nP,1]
  ests.a.se <- AAmethod.out$estimates[1:nP,2]
  m.est <- - AAmethod.out$estimates['exponent',1] # renegate
  m.est.se <- AAmethod.out$estimates['exponent',2]
  
  plot(ests.a ~ log(Ps), pch=19, ylim=c(min(ests.a-ests.a.se), max(ests.a+ests.a.se)), ylab='log(a)',xlab=('log(P)'), las=1)
  arrows(log(Ps), ests.a-ests.a.se, log(Ps), ests.a+ests.a.se, code=3, angle=90, length=0.1)
  abline(AAmethod.out$fit.m)
  legend('topright', legend=bquote(m==.(round(m.est,2))%+-%.(round(m.est.se,2))), inset=0.1, bty='n')
}

# ##############################################################
# ##############################################################
# # ~~~~~~~~~~~~~~~~~
# # Test on some data
# # ~~~~~~~~~~~~~~~~~
# # Katz 1985 data (from Arditi & Akcakya paper)
# Katz <- read.table(header=TRUE, text="
#                    Npredator	Nprey	Nconsumed
#                    1	16	2.14
#                    1	32	4.14
#                    1	64	4.29
#                    1	128	4.57
#                    2	16	1.29
#                    2	32	8.29
#                    2	64	8.14
#                    2	128	9.14
#                    3	16	1.71
#                    3	32	6.43
#                    3	64	7.86
#                    3	128	13.57
#                    4	16	2.29
#                    4	32	7.29
#                    4	64	8.71
#                    4	128	18.0
#                    ")

# # edwards_1961_Trichogramma-Sitotroga_2
# Edwards <- read.table(header=TRUE, text="
#                       Npredator	Nprey	Nconsumed
#                       1	18	1
#                       1	32	1
#                       1	72	2
#                       1	128	6
#                       1	200	2
#                       10	32	4
#                       10	72	20
#                       10	128	15
#                       10	200	11
#                       20	32	8
#                       20	72	22
#                       20	128	29
#                       20	200	37
#                       60	32	5
#                       60	72	12
#                       60	128	14
#                       60	200	30
#                       ")
# 
# 
# dropboxdir <- switch(
#   Sys.getenv("LOGNAME"),
#   stouffer = '../../../dropbox_data/Data',
#   marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
# )
# source("./Dataset_Code/Elliot_2005_Instar5.R")
# Elliot<-d

# # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat <- Katz
# dat$Time<-1
# # Note: Should be bootstrapped, but just do one draw for test
# dat$Nconsumed<-rbinom(nrow(dat),size=dat$Nprey,prob=dat$Nconsumed/dat$Nprey)
# okay4AAmethod(dat)
# katz.out <- AAmethod(dat, replacement=FALSE)
# plot.AAmethod(katz.out)
# 
# #  ~~~~~~~~~~~~~
# 
# dat <- Edwards
# okay4AAmethod(dat)
# edwards.out <- AAmethod(dat, replacement=FALSE)
# plot.AAmethod(edwards.out)

# #  ~~~~~~~~~~~~~

# dat <- Elliot
# okay4AAmethod(dat)
# elliot.out <- AAmethod(dat, replacement=TRUE)
# plot.AAmethod(elliot.out)
###################################################################
###################################################################
###################################################################


