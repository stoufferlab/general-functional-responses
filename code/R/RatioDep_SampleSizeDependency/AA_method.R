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
    Neaten <- N0 - (1 / (a * h)) * lamW::lambertW0( (a * h * N0) * exp(-a * (T - h * N0)))
  }
  
  if(expttype=="replacement"){
    numer <- (a * N0)
    denom <- (1 + a * h * N0)
    Neaten <- (numer / denom) * P * T
  }
  
  return(Neaten)
}


# New negative log likelihood
holling2.NLL = function(
  attack,
  handling,
  initial,
  killed,
  predators,
  expttype,
  Pminus1,
  time=NULL
){
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
dat <- read.table(header=TRUE, text="
                P	N	Eaten
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

MM<-model.matrix(~0+as.factor(dat$P))


a<-as.vector(1:ncol(MM))

x%*%a



AA.method <- function(N,P,Neaten)
  Plevels <- sort(unique(P))
