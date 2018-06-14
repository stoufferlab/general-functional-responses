
#############################################
# ratio-dependent-like functional responses
#############################################

# all models below take forms like:
# F = a N / P
# F = a N / P^m
# F = a N / (P + a h N)
# F = a N / (P^m + a h N)

# libraries required to use these functions
library(lamW)

# predicted number of species consumed given parameters of a ratio-dependent-like functional response
ratio.dependent.pred = function(N0, a, h, m, P, T, expttype=c("integrated","replacement")){
	expttype <- match.arg(expttype)

	if(expttype=="integrated"){
		if(h == 0){
			N <- N0 * (1 - exp(-a * T * P ^ (1 - m)))
		}else{
			Q <- P ^ m
			N <- N0 - (Q / (a * h)) * lamW::lambertW0(((a * h * N0)/ Q) * exp(- (a / Q) * (P * T - h * N0)))
		}
	}

	if(expttype=="replacement"){
		numer <- (a * N0)
		denom <- (P ^ m + a * h * N0)
		N <- (numer / denom) * P * T
	}
	
	return(N)
}

# negative log likelihood for data predicted by a ratio-dependent-like functional response
ratio.dependent.NLL = function(attack, handling, exponent, initial, killed, predators, expttype, time=NULL){
	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	attack <- exp(attack)
	handling <- exp(handling)
	exponent <- exp(exponent)

	if(is.null(time)){
		time <- 1
	}
	
	# expected number consumed
	Nconsumed <- ratio.dependent.pred(N0=initial, a=attack, h=handling, m=exponent, P=predators, T=time, expttype=expttype)

	# DEBUG if the parameters are not biologically plausible, neither should be the likelihood
	if(any(Nconsumed < 0 | is.nan(Nconsumed))){
		return(Inf)
	}

	# negative log likelihood based on proportion consumed (no replacement)
	if(expttype=="integrated"){
		nll <- -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE))
	}

	# negative log likelihood based on total number consumed (replacement)
	if(expttype=="replacement"){
		nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
	}

	return(nll)
}

# convenience function to fit a hassell-varley functional response
hv.NLL = function(params, initial, killed, predators, time, expttype){
	attack <- params[1]
	handling <- log(1E-7)
	exponent <- params[2]
	
	nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
	return(nll)
}

# convenience function to fit an arditi-ginzburg functional response
ag.NLL = function(params, initial, killed, predators, time, expttype){
	attack <- params[1]
	handling <- params[2]
	exponent <- log(1)
	
	nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
	return(nll)
}

# convenience function to fit an arditi-akcakaya functional response
aa.NLL = function(params, initial, killed, predators, time, expttype){
	attack <- params[1]
	handling <- params[2]
	exponent <- params[3]
	
	nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
	return(nll)
}
