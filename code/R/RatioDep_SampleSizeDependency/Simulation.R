# Purpose: Simulation(s) to assess the dependence of the interference-rate 
# estimator m.hat of the Arditi-Akcakaya model on treatment replication 
# (sample-size) and predator-density, and compare it to c.hat for the 
# Beddington-DeAngelis model.

# Motivation:  [copy from other script]

# Setup:
# Assume the 'true' functional response is either
    # a * N / (qr.bar * P^m + a * h * N)
    #   or
    # a * N / (1 + a * h * N + c * qr.bar * P)

# where qr.bar is the predator equivalency of an average individual predator
# averaged across treatment replicates

# but estimate each using
    # a * N / (P^m.hat + a * h * N)
    #   or
    # a * N / (1 + a * h * N + c.hat * P)

# ############################################################################## 
rm(list = ls())

# ~~~~~~~~~~~~~
# Load packages
# ~~~~~~~~~~~~~
library(sn) # For skew normal
library(tidyr) # For gather() function
library(dplyr) # For starts_with() function
library(ggplot2)
library(bbmle)
library(lamW)

# ~~~~~~~~~~~~~~~~~~~~~~
# Load fitting functions
# ~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~
# Parameter set up
# ~~~~~~~~~~~~~~~~
# Fixed functional response parameters:
  a = 0.05 #  attack rate, approximates estimate from XXXXX
  h = 0.2 # handling time, approximates estimate from XXXXX
  T = 1 # duration of experiment
  m = 1 # "interference rate" (Arditi-Akcakaya model), set to 1 for simplify comparison to m-hat
  c = 1 # "interference rate" (Beddington-DeAngelis model), set to 1 for simplify comparison to m-hat
  N = 20*(1:10)^2 # prey abundances, approximately conforms to XXXXXXXXX

# Simulation dependent variables:
  # R - number of replicates per treatment
  # P.s - predator abundances (a list, with each element corresponding to a vector of predator densities)
  # qi - the predator equivalency of the ith predator individual
  # qr - the average predator equivalency of the predator individuals in replicate r
  # qr.bar - the mean predator equivalency of an average individual predator averaged across treatment replicates
  # x.bar - location parameter of the skew normal (= qr.bar only when skew parameter alpha = 0 
  #         and hence must be calculated from given qr.bar, alpha and omega)
  
# Simulation independent variables to vary over:
  Itns = 5 # Number of iterations for each combination of R, P and case
  R.s = c(1:10,15,20,25,30,40,50,75,100,150,200,250) # Levels of replication
  P.n = 4 # Number of predator density levels
  P.b = 2:7# Log-base for predator density series
  qr.s = c(1,0.5) # desired mean of skew normal distribution 
                  # (will get converted to location parameter of skew normal distribution)
  sigma2 = c(0,0.01) # desired variance of skew normal distribution 
                  # (will get converted to scale paremeter of skew normal distribution)
  alpha = c(0,-5,5) # skew parameter of skew normal distribution (right skewed when alpha < 0)
  
# ##############################################################################  
# ~~~~~~~~~~~~~~~~~~~~~
# Convenience functions
# ~~~~~~~~~~~~~~~~~~~~~
# See Wikipedia for first two:
# To keep the variance of skew normal at var = sigma^2, solve for required omega given delta and desired sigma.
omega.sn <- function(sigma2=1, alpha=0){
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sqrt(sigma2 / (1 - ( 2 * delta^2) / pi))
  return(omega)
}
  
# To keep mean of skew normal at gr.bar, solve for required x.bar given omega and alpha and desired gr.bar.
x.bar.sn <- function(qr.bar=1, omega=0, alpha=0){
  delta <- alpha / sqrt(1 + alpha^2)
  x.bar <- qr.bar - omega * delta * (sqrt(2 / pi))
  return(x.bar)
}

# Artificially truncate the skew normal to lie within min and max.  Note that moments will change!
rtsn <- function(n=1, xi=1, omega=0, alpha=0, min=-Inf, max=Inf){
  x <- rsn(n, xi, omega, alpha)
  repl <- which(x < min | x > max)
  if(length(repl)>0){  x[repl] <- rtsn(length(repl), xi, omega, alpha, min, max) }
  return(x)
}
# z.n=100000; z.qr.bar=1; z.sigma2=0.04; z.max=2*z.qr.bar;
# z.alpha=0
#   z <- rtsn(n=z.n, xi=x.bar.sn(z.qr.bar, omega.sn(z.sigma2, z.alpha), z.alpha), omega=omega.sn(z.sigma2, z.alpha),
#               alpha=z.alpha, min=0,max=z.max); hist(z); mean(z); sd(z)
# z.alpha=5 # skew left; tail to the right
#   z <- rtsn(n=z.n, xi=x.bar.sn(z.qr.bar, omega.sn(z.sigma2, z.alpha), z.alpha), omega=omega.sn(z.sigma2, z.alpha),
#               alpha=z.alpha, min=0,max=z.max); hist(z); mean(z); sd(z)
# z.alpha=-5 # skew right; tail to the left
#   z <- rtsn(n=z.n, xi=x.bar.sn(z.qr.bar, omega.sn(z.sigma2, z.alpha), z.alpha), omega=omega.sn(z.sigma2, z.alpha),
#             alpha=z.alpha, min=0,max=z.max); hist(z); mean(z); sd(z)

# ~~~~~~~~~~~
# Logarithmic sequences (for creating predator abundance sequences of different log-base)
# (Start with P=2 so as not to run into problem of partial indivdiual when P=1)
# as.integer(round(...)) is needed because integers are needed to denote sample size (n) and R won't necessarily show that you don't have an integer, causing sample size to be rounded down!!!
log.seq<-function(base=2,length=5,min=2){as.integer(round(base^(seq(from=log(min,base),length=length))))}

# ~~~~~~~~~~~
# Implement a version of the Arditi-Akcakaya method of estimating 'm'

# ~~~~~~~~~~~~~~~~~
# Fitting functions
# ~~~~~~~~~~~~~~~~~
# library(RCurl)
# 
# script <- getURL("https://raw.githubusercontent.com/stoufferlab/general-functional-responses/master/code/R/Single_Predator_Single_Prey/fit_holling_like.R?token=AA5Y2dmhno4Umy6A0Cbu_nhxWt0hMXVFks5bKuMrwA%3D%3D", ssl.verifypeer = FALSE)
# 
# eval(parse(text = script))

# Grabbed and updated from fit_ratio_dependent.R to allow for partial predators (q parameter added)
ratio.dependent.pred = function(N0, a, h, m, q, P, T, expttype=c("integrated","replacement")){
  expttype <- match.arg(expttype)
  
  if(expttype=="integrated"){
    if(h == 0){
      N <- N0 * (1 - exp(-a * T * (q * P) ^ (1 - m)))
    }else{
      Q <- (q * P) ^ m
      N <- N0 - (Q / (a * h)) * lamW::lambertW0(((a * h * N0)/ Q) * exp(- (a / Q) * (P * T - h * N0)))
    }
  }
  
  if(expttype=="replacement"){
    numer <- (a * N0)
    denom <- ((q * P) ^ m + a * h * N0)
    N <- (numer / denom) * P * T
  }
  
  return(N)
}

ratio.dependent.NLL = function(attack, handling, exponent, partialP, initial, killed, predators, expttype, time=NULL){
  # we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
  attack <- exp(attack)
  handling <- exp(handling)
  exponent <- exp(exponent)
  partialP <- exp(partialP)

  if(is.null(time)){
    time <- 1
  }
  
  # expected number consumed
  Nconsumed <- ratio.dependent.pred(N0=initial, a=attack, h=handling, m=exponent, q=partialP, P=predators, T=time, expttype=expttype)
  
  # negative log likelihood based on proportion consumed (no replacement)
  # DEBUG: consider whether binomial or poisson are interchangeable
  if(expttype=="integrated"){
    nll <- -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE))
    # nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
  }
  
  # negative log likelihood based on total number consumed (replacement)
  if(expttype=="replacement"){
    nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
  }
  
  return(nll)
}


aa.NLL = function(params, initial, killed, predators, time, expttype){
  attack <- params[1]
  handling <- params[2]
  exponent <- params[3]
  
  nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, partialP=log(1), initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
  return(nll)
}

qaa.NLL = function(params, initial, killed, predators, time, expttype){
  attack <- params[1]
  handling <- params[2]
  exponent <- params[3]
  partialP <- params[4]
  
  nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, partialP=partialP, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
  return(nll)
}

# ############################################################################## 
# ~~~~~~~~~~~~~~~~
# Additional setup
# ~~~~~~~~~~~~~~~~
# Skew normal parameter combinations to sample from
SN <- expand.grid(qr.bar=qr.s, sigma2=sigma2, alpha=alpha)
  # Remove unneccessary combinations
  rem <- c(which(SN$sigma2==0 & SN$alpha!=0), which(SN$qr.bar!=1 & SN$sigma2!=0), which(SN$qr.bar!=1 & SN$alpha!=0))
  SN <- SN[-rem,]
  SN$omega <- omega.sn(SN$sigma2,SN$alpha)
  SN$x.bar <- x.bar.sn(SN$qr.bar,SN$omega,SN$alpha)
  SN.cases <- SN
  SN.cases$case<-1:nrow(SN)
  SN.cases$caseName <- c('1.Full.NoVar','2.Partial.NoVar','3.Full.Var','4.Full.SkewRightTailLeft','5.Full.SkewLeftTailRight')

# Create sets of predator densities
P.s <- mapply(log.seq,base=P.b,length=P.n)
  P.s <- matrix(unlist(P.s), ncol=length(P.b))
  colnames(P.s) <- P.s[nrow(P.s),] # name series according to their max predator density

# Alternative 1
P.s <- mapply(seq, from=seq(2,10,by=2), length=P.n, by=seq(2,10,by=2))
   colnames(P.s) <- P.s[nrow(P.s),] # name series according to their max predator density

# Alternative 2
P.s <- mapply(seq, from=2, length=P.n, by=seq(2,6))
colnames(P.s) <- P.s[nrow(P.s),] # name series according to their max predator density
   
# ##############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create data sets and fit models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out<-dim(0)
for(s in 1:nrow(SN)){
  for(p in 1:ncol(P.s)){
    P <- P.s[,p]
    NP <- expand.grid(P=P,N=N)
      for(R in R.s){
        for(i in 1:Itns){
        # mapply(...) -- draw P times from the truncated skew normal (distribution of qi's) for each predator level
        # lapply(...,sum) -- P.eff, the effective number of predators (mean qi across predators times the number of predators) within each predator abundance level
        # replicate(R,...) -- repeat above R times
        
        P.eff <- replicate(R, lapply(mapply(rtsn, n=round(NP$P), xi=SN$x.bar[s], omega=SN$omega[s], alpha=SN$alpha[s], min=0, max=2*SN$x.bar[s]), sum))
        P.eff <- matrix(as.numeric(P.eff), ncol=R)
        
      # Reshape data for model-fitting
        # Arditi-Akcakaya feeding rates
        Fr.AA <- apply(P.eff,2,function(x){a * NP$N * NP$P / (x^m + a * h * NP$N)})
        if(R==1){ Fr.AA <- data.frame(NP, Frate=Fr.AA) }
        if(R!=1){ Fr.AA <- gather(data.frame(NP, Fr.AA), Rep, Frate, starts_with('X')) }
        Fr.AA$Neaten <- rpois(nrow(Fr.AA),Fr.AA$Frate)
        # Fr.AA <- Fr.AA[ do.call(order, Fr.AA), ]
        # ggplot(Fr.AA,aes(N,Frate/P)) + geom_point(aes(colour=factor(P)))
        # ggplot(Fr.AA,aes(N,Neaten/P)) + geom_point(aes(colour=factor(P)))

        # Beddington-DeAngelis feeding rates
        Fr.BD <- apply(P.eff,2,function(x){a * NP$N * NP$P / (1 + a * h * NP$N + c * x)})
        if(R==1){ Fr.BD <- data.frame(NP, Frate=Fr.BD) }
        if(R!=1){ Fr.BD <- gather(data.frame(NP,Fr.BD), Rep, Frate, starts_with('X')) }
        Fr.BD$Neaten <- rpois(nrow(Fr.BD),Fr.BD$Frate)
        # Fr.BD <- Fr.BD[ do.call(order, Fr.BD), ]
        # ggplot(Fr.BD,aes(N,Frate/P)) + geom_point(aes(colour=factor(P)))

      # Fit models
        # Arditi-Akcakaya
        fit.AA.opt <- nloptr::nloptr(
          x0=c( log(a), log(h), log(m) ),
          eval_f=aa.NLL,
          opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
          initial=Fr.AA$N,
          killed=Fr.AA$Neaten,
          predators=Fr.AA$P,
          time=1,
          expttype='replacement'
        )
        st.AA<-fit.AA.opt$solution; names(st.AA)<-c('attack','handling','exponent')
        
        fit.AA <- bbmle::mle2(
          ratio.dependent.NLL,
          start=list(attack=st.AA[1], handling=st.AA[2], exponent=st.AA[3]),
          fixed=list(partialP=log(1)),
          data=list(initial=Fr.AA$N, killed=Fr.AA$Neaten, predators=Fr.AA$P, 
                    time=1, expttype='replacement', params=NULL))
        # print(summary(fit.AA))
        est.AA <- coef(fit.AA); names(est.AA) <- paste0('AA.',names(est.AA))
        
        
        # Arditi-Akcakaya allowing for partial predators
        fit.qAA.opt <- nloptr::nloptr(
          x0=c( log(a), log(h), log(m), log(1) ),
          eval_f=qaa.NLL,
          opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
          initial=Fr.AA$N,
          killed=Fr.AA$Neaten,
          predators=Fr.AA$P,
          time=1,
          expttype='replacement'
        )
        st.qAA<-fit.qAA.opt$solution; names(st.qAA)<-c('attack','handling','exponent','partialP')
        
        fit.qAA <- bbmle::mle2(
          ratio.dependent.NLL,
          start=list(attack=st.qAA[1], handling=st.qAA[2], exponent=st.qAA[3], partialP=st.qAA[4]),
          data=list(initial=Fr.AA$N, killed=Fr.AA$Neaten, predators=Fr.AA$P, 
                    time=1, expttype='replacement', params=NULL))
        est.qAA <- coef(fit.qAA); names(est.qAA) <- paste0('qAA.',names(est.qAA))
        
        # Beddington-DeAngelis
        # fit.BD

        
        # Combine output
        out <- rbind(out, 
              t(as.matrix(c(a=a, h=h, m=m, c=c, SN[s,], mxP=max(NP$P), R=R, Itn=i, est.AA, est.qAA))))
        
 print(paste0('Case ',s, '; Reps ',R, '; maxP ',max(NP$P), '; Iteration ',i,' of ',Itns))
 
}}}}
dat <- data.frame(matrix(unlist(out),ncol=ncol(out))); colnames(dat) <- colnames(out)

# Back-transform parameters
dat[,-c(1:12)] <- apply(dat[,-c(1:12)],2,exp)
dat <- merge(SN.cases,dat) # Merge case names in

# Export results
save(dat,file='Simulation_dat.Rdata')

# ##############################################################################
# ~~~~~~~~~~~~
# Plot results
# ~~~~~~~~~~~~
load(file='Simulation_dat.Rdata')

# Rearrange models to enable use of facets
dat.AA <- dat[,c(1:14,15:18)]
  colnames(dat.AA)[15:18] <- sub('AA.','',colnames(dat.AA)[15:18])
  dat.AA$model <- 'AA'
dat.qAA <- dat[,c(1:14,19:22)]    
  colnames(dat.qAA)[15:18] <- sub('qAA.','',colnames(dat.qAA)[15:18])
  dat.qAA$model <- 'qAA'
dat2 <- rbind(dat.AA,dat.qAA)
               

ggplot(dat2,aes(R,exponent, colour=factor(mxP))) + geom_point() + geom_smooth(se=FALSE) + facet_grid(model~caseName) + coord_trans(x="log")
ggplot(dat2,aes(R,attack, colour=factor(mxP))) + geom_point() + geom_smooth(se=FALSE) + facet_grid(model~caseName) + coord_trans(x="log")
ggplot(dat2,aes(R,handling, colour=factor(mxP))) + geom_point() + geom_smooth(se=FALSE) + facet_grid(model~caseName) + coord_trans(x="log")

ggplot(dat2,aes(attack,exponent)) + geom_point(aes(colour=factor(mxP))) + facet_grid(model~caseName)
ggplot(dat2,aes(handling,exponent)) + geom_point(aes(colour=factor(mxP))) + facet_grid(model~caseName)
ggplot(dat2,aes(attack,handling)) + geom_point(aes(colour=factor(mxP))) + facet_grid(model~caseName)

