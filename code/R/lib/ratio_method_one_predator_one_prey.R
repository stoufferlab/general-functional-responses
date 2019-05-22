#############################################
# ratio-dependent-like functional responses
# modified from 'holling-like functional responses'
#############################################
# For non-replacement datasets, we have two options : 
# (1) solve using lambertsW (or solve transcendental eqn directly if needed)
# (2) integrate
#############################################
# libraries required
library(bbmle)
library(nloptr)
library(lamW)
library(odeintr)
library(emdbook) # for beta-binomial 

sp <- list.files("../../..", "set_params.R", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
source(sp)
#############################################

# For integration method, define the ode in C++ format
# ratio.like.1pred.1prey.sys = '
#   // ratio-dependent-family of functional responses for one predator one prey
#   dxdt[0] = -P * (a * x[0]) / (pow(P , m) + a * h * x[0]);
#   
#   // consumption rate cannot be positive
#   if(dxdt[0] > 0) dxdt[0] = 0;
# '
# 
# # compile the above C++ code into something we can run in R
# odeintr::compile_sys(
#   "ratio_1pred_1prey",
#   ratio.like.1pred.1prey.sys,
#   pars = c("a", "h", "m", "P") #,
#   # method = "bsd"
# )

# predicted number of species consumed given parameters of a ratio-dependent family functional response
ratio.like.1pred.1prey = function(N0, a, h, m, P, T, 
                                  replacement, 
                                  integrate=FALSE,
                                  overrideTranscendental=FALSE){

	# in a world with replacement everything is hunky dory
	if(replacement){
	  numer <- (a * N0)
	  denom <- (P ^ m + a * h * N0)
	  N <- (numer / denom) * P * T
	  N <- pmax(0,N)
	  return(N)
	}

	# without replacement
	if(!replacement){
	  if(h==0){ # For Hassell-Varley ("Type 1") things are simple:
	    N <- N0 * (1 - exp(-a * T * P ^ (1 - m)))
		}else{# For all other models...
		  if(integrate){  # solve by direct integration
		    N <- numeric(length(N0))
		    for(i in seq.int(length(N0))){
		      
		      # set parameters within ode solver
		      ratio_1pred_1prey_set_params(a=a, h=h, m=m, P=P[i])
		      
		      # calculate the final number of prey integrating the ode
		      Nfinal <- ratio_1pred_1prey(N0[i], T[i], T[i]/100.)
		      
		      # we only need the last row since this is the final "abundance"
		      Nfinal <- as.numeric(Nfinal[nrow(Nfinal),2])
		      
		      # the number consumed is the difference between what we started with and what is left
		      N[i] <- N0[i] - Nfinal
		    }
		  } else{  # solve using lambertsW (or trancendental equation)
    		  Q <- P ^ m
    		  N <- N0 - (Q / (a * h)) * lamW::lambertW0(((a * h * N0) / Q) * exp(- (a / Q) * (P * T - h * N0)))
    			
    			# sometimes the argument in the exponential passed to lambertW0 causes it to blow up
    		  if(!overrideTranscendental){
      			if(any(is.infinite(N))){
      				# the explicit result of the analytical integration without solving for N implictly
      				ffff <- function(N, N0, P, T, a, h, Q){
      					dN <- Q * log((N0 - N) / N0) - a * h * N
      					dt <- - a * P * T
      					dN - dt
      				}
      				# sometimes the time argument is a constant and not a vector
      				if(length(T)==1){
      					T <- rep(T, length(N0))
      				}
      				# check which predictions are non-sensical
      				for(i in 1:length(N0)){
      					if(is.infinite(N[i])){
      						# we need to solve the transcendental equation directly
      						nn <- uniroot(ffff, lower=0, upper=N0[i], N0=N0[i], P=P[i], T=T[i], a=a, h=h, Q=Q[i])
      						N[i] <- nn$root
      					}
      				}
      			}
    		  }
		  }
		}
		return(N)
	}

	stop()
}

# negative log likelihood for ratio-like models given parameters and requisite data
ratio.like.1pred.1prey.NLL = function(params, 
                                      modeltype, 
                                      initial, 
                                      killed, 
                                      predators, 
                                      replacement, 
                                      time=NULL){
  
	set_params(params, modeltype)


	# if no times are specified then normalize to time=1
	if(is.null(time)){
		time <- 1
	}

	# expected number consumed given data and parameters
	Nconsumed <- ratio.like.1pred.1prey(N0=initial,
	                                    a=attack,
	                                    h=handling,
	                                    m=exponent,
	                                    P=predators,
	                                    T=time,
	                                    replacement=replacement)
  
  # reduce to unique data rows to speed up. There's probably an even faster way, but...
  # d.ori <- data.frame(predators, initial, time)
  # d.ori$id  <- 1:nrow(d.ori) # needed to reorder predictions after merge
  # d.uniq <- unique(d.ori)
  # Nconsumed.uniq <- ratio.like.1pred.1prey(N0 = d.uniq$initial,
  #                                          a = attack,
  #                                          h = handling,
  #                                          m = exponent,
  #                                          P = d.uniq$predators,
  #                                          T = d.uniq$time,
  #                                          replacement = replacement)
  # d.pred <- merge(d.ori, cbind(d.uniq, Nconsumed.uniq), all.x=TRUE)
  # Nconsumed <- d.pred$Nconsumed.uniq[order(d.pred$id)]

	# if the parameters are not biologically plausible, neither should be the likelihood
	if(any(Nconsumed <= 0) | any(is.nan(Nconsumed))){
		nll <- Inf
		return(nll)
	}else{
		# negative log likelihood based on proportion consumed (no replacement)
		if(!replacement){
			# nll <-  -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE))
			nll <- -sum(dbetabinom(killed, prob=Nconsumed/initial, size=initial, theta=theta, log=TRUE))
			if(is.nan(nll)){
			  nll <- Inf
			}
			return(nll)
		}

		# negative log likelihood based on total number consumed (replacement)
		if(replacement){
			# nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
			nll <- -sum(dnbinom(killed, mu=Nconsumed, size=theta, log=TRUE))
			return(nll)
		}
	}

	stop()
}

# needed by mle2 to pass named parameters in the right order
parnames(ratio.like.1pred.1prey.NLL) <- c(
  'theta',
	'attack',
	'handling',
	'exponent'
)

# given data (d), study info (s), and modeltype (e.g., "Ratio I"), fit functional response data
fit.ratio.like <- function(d, s, 
                           modeltype, 
                           nloptr.control=list(), 
                           mle2.control=list(), 
                           ...){
  
	# estimate starting value from the data using linear regression
  start <- list(
    theta = log(100),
    attack = log(coef(lm(d$Nconsumed~0+I(d$Npredator * d$Nprey))))
  )
	
	# fit Ratio ("Type I") via MLE with above starting parameter value
	ratio.via.sbplx <- nloptr::sbplx(
		x0 = unlist(start),
		fn = ratio.like.1pred.1prey.NLL,
		modeltype="Ratio",
		initial=d$Nprey,
		killed=d$Nconsumed,
		predators=d$Npredator,
		time=d$Time,
		replacement=s$replacement,
		control = nloptr.control,
		...
	)

	# refit with mle2 since this also estimates the covariance matrix for the parameters
	mle2.start <- as.list(ratio.via.sbplx$par)
	names(mle2.start) <- names(start)
	
	ratio.via.mle2 <- bbmle::mle2(
		minuslogl = ratio.like.1pred.1prey.NLL,
		start = mle2.start,
		data = list(
			modeltype="Ratio",
			initial=d$Nprey,
			killed=d$Nconsumed,
			predators=d$Npredator,
			time=d$Time,
			replacement=s$replacement
		),
		vecpar = TRUE,
		control = mle2.control,
		...
	)
	
	if(modeltype == "Ratio"){
	  return(ratio.via.mle2)
	}
	# code to fit subsequent models
	else{
	  # fit Arditi-Ginzburg (Ratio type II) first
	  start <- list(
	    theta = log(100),
	    attack = coef(ratio.via.mle2)["attack"],
	    handling = log(1)
	  )
	  
	  # fit the more complex model with sbplx first
	  ag.via.sbplx <- nloptr::sbplx(
	    x0 = unlist(start),
	    fn = ratio.like.1pred.1prey.NLL,
	    modeltype = "Arditi.Ginzburg",
	    initial = d$Nprey,
	    killed = d$Nconsumed,
	    predators = d$Npredator,
	    time = d$Time,
	    replacement = s$replacement,
	    control = nloptr.control,
	    ...
	  )
	  
	  mle2.start <- as.list(ag.via.sbplx$par)
	  names(mle2.start) <- names(start)
	  
	  # refit with mle2 since this also estimates the covariance matrix for the parameters
	  ag.via.mle2 <- bbmle::mle2(
	    minuslogl = ratio.like.1pred.1prey.NLL,
	    start = mle2.start,
	    data = list(
	      modeltype = "Arditi.Ginzburg",
	      initial = d$Nprey,
	      killed = d$Nconsumed,
	      predators = d$Npredator,
	      time = d$Time,
	      replacement = s$replacement
	    ),
	    vecpar = TRUE,
	    control = mle2.control
	  )
	  
	  if(modeltype == "Arditi.Ginzburg"){
	    return(ag.via.mle2)
	  }else{
	    if(modeltype == "Hassell.Varley"){
	      start <- list(
	        theta = log(100),
	        attack = coef(ratio.via.mle2)["attack"],
	        # handling = 0,
	        exponent = log(1)
	      )
	    }
	    
	    if(modeltype == "Arditi.Akcakaya"){
	      start <- list(
	        theta = log(100),
	        attack = coef(ag.via.mle2)["attack"],
	        handling = log(1),
	        exponent = log(1)
	      )
	    }
	    

	    # fit the more complex model with sbplx first
	    fit.via.sbplx <- nloptr::sbplx(
	      x0 = unlist(start),
	      fn = ratio.like.1pred.1prey.NLL,
	      modeltype = modeltype,
	      initial = d$Nprey,
	      killed = d$Nconsumed,
	      predators = d$Npredator,
	      time = d$Time,
	      replacement = s$replacement,
	      control = nloptr.control,
	      ...
	    )
	    
	    # convert nloptr estimation to list of starting values
	    mle2.start <- as.list(fit.via.sbplx$par)
	    names(mle2.start) <- names(start)
	    
	    # fit with mle2 since this provides other convenience estimates
	    fit.via.mle2 <- bbmle::mle2(
	      minuslogl = ratio.like.1pred.1prey.NLL,
	      start = mle2.start,
	      data = list(
	        modeltype = modeltype,
	        initial = d$Nprey,
	        killed = d$Nconsumed,
	        predators = d$Npredator,
	        time = d$Time,
	        replacement = s$replacement
	      ),
	      vecpar = TRUE,
	      control = mle2.control
	    )
	    
	    
	    # convert mle2 estimation to list of starting values
	    mle2.start <- as.list(fit.via.mle2@coef)
	    names(mle2.start) <- names(start)
	    
	    # apparently this helps optimize over complex likelihood surfaces and get SEs when they weren't there otherwise...
	    mle2.control$parscale <- abs(fit.via.mle2@coef)
	    
	    # refit with mle2 using parscale to help get an appropriate covariance matrix for the parameters
	    refit.via.mle2 <- bbmle::mle2(
	      minuslogl = ratio.like.1pred.1prey.NLL,
	      start = mle2.start,
	      data = list(
	        modeltype = modeltype,
	        initial = d$Nprey,
	        killed = d$Nconsumed,
	        predators = d$Npredator,
	        time = d$Time,
	        replacement = s$replacement
	      ),
	      vecpar = TRUE,
	      control = mle2.control				
	    )
	    
	    return(refit.via.mle2)
	    
	    
	  }
	}
}
