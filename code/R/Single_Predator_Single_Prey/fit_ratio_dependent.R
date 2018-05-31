
# libraries used below
library(bbmle)
library(lamW)
library(nloptr)

###############################
# ratio-dependent-like functional responses
###############################

# F = a N / P
# F = a N / P^m
# F = a N / (P + a h N)
# F = a N / (P^m + a h N)

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

ratio.dependent.NLL = function(attack, handling, exponent, initial, killed, predators, time, expttype){
	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	attack <- exp(attack)
	handling <- exp(handling)
	exponent <- exp(exponent)
	
	# expected number consumed
	Nconsumed <- ratio.dependent.pred(N0=initial, a=attack, h=handling, m=exponent, P=predators, T=time, expttype=expttype)

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

hv.NLL = function(params, initial, killed, predators, time, expttype){
	attack <- params[1]
	handling <- -20
	exponent <- params[2]
	
	nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
	return(nll)
}

ag.NLL = function(params, initial, killed, predators, time, expttype){
	attack <- params[1]
	handling <- params[2]
	exponent <- 0
	
	nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
	return(nll)
}

aa.NLL = function(params, initial, killed, predators, time, expttype){
	attack <- params[1]
	handling <- params[2]
	exponent <- params[3]
	
	nll <- ratio.dependent.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
	return(nll)
}

#####################################
# series of functional response fits
#####################################

ffr.hv.nloptr <- nloptr::nloptr(
	x0=x0.rd,
	eval_f=hv.NLL,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
)

ffr.hv <- bbmle::mle2(
	ratio.dependent.NLL,
	start=list(attack=ffr.hv.nloptr$solution[1], exponent=ffr.hv.nloptr$solution[2]),
	fixed=list(handling=-20),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, params=NULL),
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.ag.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hv)[1],-10),
	eval_f=ag.NLL,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
)

ffr.ag <- bbmle::mle2(
	ratio.dependent.NLL,
	start=list(attack=ffr.ag.nloptr$solution[1], handling=ffr.ag.nloptr$solution[2]),
	fixed=list(exponent=0),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, params=NULL),
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.aa.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.ag)[1],coef(ffr.ag)[2],0),
	eval_f=aa.NLL,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
)

ffr.aa <- bbmle::mle2(
	ratio.dependent.NLL,
	start=list(attack=ffr.aa.nloptr$solution[1], handling=ffr.aa.nloptr$solution[2], exponent=ffr.aa.nloptr$solution[3]),
	# fixed=list(exponent=-20),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, params=NULL),
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)