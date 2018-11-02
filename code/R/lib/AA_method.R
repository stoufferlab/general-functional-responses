# Implement 'Method II' of Arditi & Akcakaya 1990.
# The method fits a type II functional response to each predator level assuming a common handling time (i.e. it estimates an attack rate for each level), then regresses the log-transformed attack rates on log-transformed predator levels to estimate 'm', weighting each predator level by the inverse of its attack rate's uncertainty (variance).
#######################################################################
require(lamW)
require(bbmle)
require(nloptr)
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Holling Type II functional response
# (Simplified from other likelihood scripts)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predicted number of prey consumed
holling2 = function(N0, a, h, P, T, replacement){
  
  if(!replacement){
    Neaten <- N0 - (1 / (a * h)) * lamW::lambertW0( (a * h * N0) * exp(-a * (P *T - h * N0)))
  }
  
  if(replacement){
    numer <- (a * N0)
    denom <- (1 + a * h * N0)
    Neaten <- (numer / denom) * P * T
  }
  
  return(Neaten)
}

# ~~~~~~~~~~~~~~~~~~~~~~~
# Negative log likelihood
# ~~~~~~~~~~~~~~~~~~~~~~~
AAM.NLL = function(
  params,
  initial,
  killed,
  predators,
  replacement,
  time=NULL
){
  
  # assuming P-specific attack rates first and handling time last...
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
  
  # expected number consumed given data and parameters
  Nconsumed <- holling2(N0=initial, a=attack, h=handling, P=predators, T=time, replacement=replacement)
  
  # DEBUG if the parameters are not biologically plausible, neither should be the likelihood
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
  
  fit.AAM.nloptr <- nloptr::nloptr(
    x0=c(runif(nP),log(1)), # random starting values
    eval_f=AAM.NLL,
    opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
    initial=d$Nprey,
    killed=d$Nconsumed,
    predators=d$Npredator,
    time=NULL,
    replacement=replacement
  )
  
  AAM.start<-fit.AAM.nloptr$solution
  names(AAM.start) <- parnames(AAM.NLL) <- c(paste0('a',1:nP),'h')
  
  # Refit with mle with nloptr starting estimates 
  # (for convenience to get errors)
  fit.AAM.mle <- bbmle::mle2(
    AAM.NLL,
    start=AAM.start,
    data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, replacement=replacement)
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Estimate 'm' as slope of log(a's) ~ log(P's)
  # using "the reciprocal of the variance of y as the weight w."
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ests <- coef(summary(fit.AAM.mle))
  ests.a <- ests[1:nP,1]
  Ps <- unique(attributes(fit.AAM.mle)$data$predators) # grab from fit to ensure Ps are in order corresponding to ests
  w <- 1/diag(vcov(fit.AAM.mle))[1:nP]   # Regn weights

  # Estimates are already log-transformed by likelihood function
  fit.AAM.lm <- lm(ests.a ~ log(Ps), weights=w)
  ests.m <- coef(summary(fit.AAM.lm))
    rownames(ests.m) <- c('a0','m')
  
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
  m.est <- AAmethod.out$estimates['m',1]
  m.est.se <- AAmethod.out$estimates['m',2]
  
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
# 
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
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat <- Katz
# # Note: Should be bootstrapped, but just do one draw for test
# dat$Nconsumed<-rbinom(nrow(dat),size=dat$Nprey,prob=dat$Nconsumed/dat$Nprey)
# 
# okay4AAmethod(dat)
# katz.out <- AAmethod(dat, replacement=FALSE)
# plot.AAmethod(katz.out)
# 
# # ~~~~~~~~~~~~~
# dat <- Edwards
# okay4AAmethod(dat)
# edwards.out <- AAmethod(dat, replacement=FALSE)
# plot.AAmethod(edwards.out)

###################################################################
###################################################################
###################################################################


