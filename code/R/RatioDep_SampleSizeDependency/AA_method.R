# Implement 'Method II' of Arditi & Akcakaya 1990

# Assumes values of P represent discrete levels, with variation in N for each P 
# (i.e. a 'manipulative experiment'), rather than a set of continuous densities 

library(lamW)

#####################################
# Holling Type II functional response
#####################################
# Predicted number of prey consumed
holling2 = function(N0, a, h, P, T, expttype=c("integrated","replacement")){
  expttype <- match.arg(expttype)
  
  if(expttype=="integrated"){
    Neaten <- N0 - (1 / (a * h)) * lamW::lambertW0( (a * h * N0) * exp(-a * (P *T - h * N0)))
  }
  
  if(expttype=="replacement"){
    numer <- (a * N0)
    denom <- (1 + a * h * N0)
    Neaten <- (numer / denom) * P * T
  }
  
  return(Neaten)
}


# Negative log likelihood for mle2
AAM.NLL = function(
  params,
  initial,
  killed,
  predators,
  expttype,
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
  Nconsumed <- holling2(N0=initial, a=attack, h=handling, P=predators, T=time, expttype=expttype)
  
  # DEBUG if the parameters are not biologically plausible, neither should be the likelihood
  if(any(Nconsumed < 0 | is.nan(Nconsumed))){
    return(Inf)
  }
  
  # negative log likelihood based on proportion consumed (no replacement)
  # DEBUG: consider whether binomial or poisson are interchangeable
  if(expttype=="integrated"){
    nll <- -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE))
  }
  
  # negative log likelihood based on total number consumed (replacement)
  if(expttype=="replacement"){
    nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
  }
  
  return(nll)
}


###############################################################################
# Grab the Katz 1985 data (from Arditi & Akcakya paper)
Katz <- read.table(header=TRUE, text="
                P	N	Neaten
                1	16	2.14
                1	32	4.14
                1	64	4.29
                1	128	4.57
                2	16	1.29
                2	32	8.29
                2	64	8.14
                2	128	9.14
                3	16	1.71
                3	32	6.43
                3	64	7.86
                3	128	13.57
                4	16	2.29
                4	32	7.29
                4	64	8.71
                4	128	18.0
                ")
# Should bootstrap the following.
dat <- Katz
dat$Neaten<-rbinom(nrow(dat),size=dat$N,prob=dat$Neaten/dat$N)


# Grab edwards_1961_Trichogramma-Sitotroga_2 since we get a high 'm' estimate for it (~1.768)
Edwards <- read.table(header=TRUE, text="
              P	N	Neaten
              1	18	1
              1	32	1
              1	72	2
              1	128	6
              1	200	2
              10	32	4
              10	72	20
              10	128	15
              10	200	11
              20	32	8
              20	72	22
              20	128	29
              20	200	37
              60	32	5
              60	72	12
              60	128	14
              60	200	30
              ")
dat <- Edwards
###############################################################################
######################
# Fit to the Katz data
######################
nP <- length(unique(dat$P))

fit.AAM.nloptr <- nloptr::nloptr(
  x0=c(runif(nP),log(1)),
  eval_f=AAM.NLL,
  opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
  initial=dat$N,
  killed=dat$Neaten,
  predators=dat$P,
  time=NULL,
  expttype="integrated"
)

AAM.start<-fit.AAM.nloptr$solution
names(AAM.start) <- parnames(AAM.NLL) <- c(paste0('a',1:nP),'h')

fit.AAM.mle <- bbmle::mle2(
  AAM.NLL,
  start=AAM.start,
  data=list(initial=dat$N, killed=dat$Neaten, predators=dat$P, expttype='integrated')
)

summary(fit.AAM.mle)

##############################################################
# Estimate 'm' as slope of log(a's) ~ log(P's)
# using "the reciprocal of the variance of y as the weight w."
##############################################################
ests <- coef(summary(fit.AAM.mle)) 
ests.a <- ests[-nrow(ests),1] # Are already log-transformed, so no need to do again
ests.a.se <- ests[-nrow(ests),2]
Ps <- unique(attributes(fit.AAM.mle)$data$predators) # grab from fit to ensure Ps are in order corresponding to ests
w <- diag(vcov(fit.AAM.mle))

# Estimate 'm' as slope
fit.AAM.lm <- lm(ests.a ~ log(Ps), weights=w[-nrow(ests)])
summary(fit.AAM.lm)

m.est <- coef(fit.AAM.lm)[2]
m.est.se  <- coef(summary(fit.AAM.lm))[2,2]



plot(ests.a ~ log(Ps), pch=19, ylim=c(min(ests.a-ests.a.se),max(ests.a+ests.a.se)),ylab='log(a)',xlab=('log(P)'))
  arrows(log(Ps), ests.a-ests.a.se, log(Ps), ests.a+ests.a.se, code=3, angle=90,length=0.1)
  abline(fit.AAM.lm)
  legend('topright',legend=bquote(m == .(round(m.est,2)) %+-%  .(round(m.est.se,2))),inset=0.1,bty='n')

##############################################################
