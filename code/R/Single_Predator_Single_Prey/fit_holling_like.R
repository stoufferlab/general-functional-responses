
# libraries used below
library(bbmle)
library(lamW)
library(nloptr)

###############################
# generalizable functional responses
###############################

# WARNING: consider the P-1 implications below!
general.pred = function(N0, a, h, c, phi_numer, phi_denom, P, T, expttype=c("integrated","replacement"), Pminus1=c(TRUE,FALSE)){
	expttype <- match.arg(expttype)
	# Pminus1 <- match.arg(Pminus1)

	# DEBUG make sure this carries through	
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

general.NLL = function(attack, handling, interference, phi_numer, phi_denom, initial, killed, predators, expttype, Pminus1, time=NULL){
	# DEBUG we should probably force some of these depending on model type
	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	attack <- exp(attack)
	handling <- exp(handling)
	interference <- exp(interference)
	phi_numer <- plogis(phi_numer)
	phi_denom <- plogis(phi_denom)

	if(is.null(time)){
		time <- 1
	}

	# expected number consumed
	Nconsumed <- general.pred(N0=initial, a=attack, h=handling, c=interference, phi_numer=phi_numer, phi_denom=phi_denom, P=predators, T=time, expttype=expttype, Pminus1=Pminus1)

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

general.NLL.params = function(params, modeltype, initial, killed, predators, time, expttype, Pminus1){
	if(modeltype == "Holling I"){
		attack <- params[1]
		handling <- log(1E-7)
		interference <- log(1E-7)
		phi_numer <- -1 * qlogis(1E-9)
		phi_denom <- -1 * qlogis(1E-9)
	}

	if(modeltype == "Holling II"){
		attack <- params[1]
		handling <- params[2]
		interference <- log(1E-7)
		phi_numer <- -1 * qlogis(1E-9)
		phi_denom <- -1 * qlogis(1E-9)
	}

	if(modeltype == "Beddington-DeAngelis"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- -1 * qlogis(1E-9)
		phi_denom <- -1 * qlogis(1E-9)
	}

	if(modeltype == "Crowley-Martin"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- -1 * qlogis(1E-9)
		phi_denom <- qlogis(1E-9)
	}

	if(modeltype == "Stouffer-Novak I"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- -1 * qlogis(1E-9)
		phi_denom <- params[4]
	}

	if(modeltype == "Stouffer-Novak II"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- params[4]
		phi_denom <- -1 * qlogis(1E-9)
	}
	
	if(modeltype == "Stouffer-Novak III"){
		attack <- params[1]
		handling <- params[2]
		interference <- params[3]
		phi_numer <- params[4]
		phi_denom <- params[5]
	}

	nll <- general.NLL(attack=attack, handling=handling, interference=interference, phi_numer=phi_numer, phi_denom=phi_denom, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype, Pminus1=Pminus1)
	return(nll)
}

#####################################
# series of functional response fits
#####################################

ffr.hollingI.nloptr <- nloptr::nloptr(
	x0=x0.hl,
	eval_f=general.NLL.params,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Holling I",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

ffr.hollingI <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.hollingI.nloptr$solution[1]),
	fixed=list(
		handling = log(1E-7),
		interference = log(1E-7),
		phi_numer = -1 * qlogis(1E-9),
		phi_denom = -1 * qlogis(1E-9)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.hollingII.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingI)[1],log(1)),
	eval_f=general.NLL.params,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Holling II",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

ffr.hollingII <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.hollingII.nloptr$solution[1], handling=ffr.hollingII.nloptr$solution[2]),
	fixed=list(
		interference = log(1E-7),
		phi_numer = -1 * qlogis(1E-9),
		phi_denom = -1 * qlogis(1E-9)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL),
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.bd.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingII)[1:2],log(1)),
	eval_f=general.NLL.params,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Beddington-DeAngelis",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

ffr.bd <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.bd.nloptr$solution[1], handling=ffr.bd.nloptr$solution[2], interference=ffr.bd.nloptr$solution[3]),
	fixed=list(
		phi_numer = -1 * qlogis(1E-9),
		phi_denom = -1 * qlogis(1E-9)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
)

ffr.cm.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingII)[1:2],log(1)),
	eval_f=general.NLL.params,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Crowley-Martin",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

ffr.cm <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.cm.nloptr$solution[1], handling=ffr.cm.nloptr$solution[2], interference=ffr.cm.nloptr$solution[3]),
	fixed=list(
		phi_numer = -1 * qlogis(1E-9),
		phi_denom = qlogis(1E-9)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL),
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.sn1.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingII)[1:2],log(1),qlogis(0.5)),
	eval_f=general.NLL.params,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Stouffer-Novak I",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

ffr.sn1 <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.sn1.nloptr$solution[1], handling=ffr.sn1.nloptr$solution[2], interference=ffr.sn1.nloptr$solution[3], phi_denom=ffr.sn1.nloptr$solution[4]),
	fixed=list(
		phi_numer = -1 * qlogis(1E-9)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
)

# DEBUG keeps failing for some datasets
# ffr.sn2.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],0,0),
# 	eval_f=sn1.NLL,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# )

# ffr.sn2 <- bbmle::mle2(
# 	general.NLL,
# 	start=list(attack=ffr.sn2.nloptr$solution[1], handling=ffr.sn2.nloptr$solution[2], interference=ffr.sn2.nloptr$solution[3], phi_numer=ffr.sn2.nloptr$solution[4]),
# 	fixed=list(phi_denom=20),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype),
# )

ffr.sn3.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingII)[1:2],log(1),qlogis(0.5),qlogis(0.5)),
	eval_f=general.NLL.params,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Stouffer-Novak III",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

ffr.sn3 <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.sn3.nloptr$solution[1], handling=ffr.sn3.nloptr$solution[2], interference=ffr.sn3.nloptr$solution[3], phi_numer=ffr.sn3.nloptr$solution[4], phi_denom=ffr.sn3.nloptr$solution[5]),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
)
