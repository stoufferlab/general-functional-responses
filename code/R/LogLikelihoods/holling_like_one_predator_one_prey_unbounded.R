
#############################################
# holling-like functional responses
#############################################

# libraries required below
library(lamW)

# predicted number of species consumed given parameters of a holling-like functional response
holling.like.1pred.1prey = function(N0, a, h, c, phi_numer, phi_denom, P, T, expttype=c("integrated","replacement"), Pminus1=c(TRUE,FALSE)){
	expttype <- match.arg(expttype)

	# if only P-1 individuals interference with predators that are doing the feeding
	if(Pminus1){
		P_interfering <- P - 1
	}else{
		P_interfering <- P
	}

	if(expttype=="integrated"){
		h <- h * (1 + (1 - phi_numer * phi_denom) * a * h * c * P_interfering)
		Q <- (1 + c * P_interfering)
		X <- (1 + (1 - phi_numer) * c * P_interfering)
		N <- N0 - (Q / (a * h)) * lamW::lambertW0(((a * h * N0)/ Q) * exp(- (a / Q) * (X * P * T - h * N0)))
	}

	if(expttype=="replacement"){
		numer <- (a * N0 * (1 + (1 - phi_numer) * c * P_interfering))
		denom <- (1 + a * h * N0 + c * P_interfering + (1 - phi_numer * phi_denom) * c * P_interfering * a * h * N0)
		N <- (numer / denom) * P * T
	}

	return(N)
}

# negative log likelihood for data predicted by a ratio-dependent-like functional response
holling.like.1pred.1prey.NLL = function(
	attack,
	handling,
	interference,
	phi_numer,
	phi_denom,
	initial,
	killed,
	predators,
	expttype,
	Pminus1,
	# phi_link=c('identity', 'logit'),
	time=NULL
){
	# phi_link <- match.arg(phi_link)
	# DEBUG we should probably force some of these here depending on model type instead of (below)

	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	attack <- exp(attack)
	handling <- exp(handling)
	interference <- exp(interference)

	# # at times we will constrain these parameters to be [0,1]
	# phi_numer <- switch(phi_link,
	# 	identity = identity(phi_numer),
	# 	logit = plogis(phi_numer)
	# )
	# phi_denom <- switch(phi_link,
	# 	identity = identity(phi_denom),
	# 	logit = plogis(phi_denom)
	# )

	# if no times are specified then normalize to time=1
	if(is.null(time)){
		time <- 1
	}

	# expected number consumed given data and parameters
	Nconsumed <- holling.like.1pred.1prey(N0=initial, a=attack, h=handling, c=interference, phi_numer=phi_numer, phi_denom=phi_denom, P=predators, T=time, expttype=expttype, Pminus1=Pminus1)

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

# convenience function for these holling-like models that will allow the use of nloptr instead of mle2
holling.like.1pred.1prey.NLL.params = function(params, modeltype, initial, killed, predators, expttype, Pminus1, time=NULL){
	if(modeltype == "Holling I"){
		attack <- params[1]
		handling <- log(1E-7)
		interference <- log(1E-7)
		phi_numer <- 1
		phi_denom <- 1
	}

	if(modeltype == "Holling II"){
		attack <- params[1]
		handling <- params[2]
		interference <- log(1E-7)
		phi_numer <- 1
		phi_denom <- 1
	}

	if(modeltype == "Beddington-DeAngelis"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- 1
		phi_denom <- 1
	}

	if(modeltype == "Crowley-Martin"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- 1
		phi_denom <- 0
	}

	if(modeltype == "Stouffer-Novak I"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- 1
		phi_denom <- params[4]
	}

	if(modeltype == "Stouffer-Novak II"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- params[4]
		phi_denom <- 1
	}
	
	if(modeltype == "Stouffer-Novak III"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- params[4]
		phi_denom <- params[5]
	}

	nll <- holling.like.1pred.1prey.NLL(attack=attack, handling=handling, interference=interference, phi_numer=phi_numer, phi_denom=phi_denom, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype, Pminus1=Pminus1)
	return(nll)
}
