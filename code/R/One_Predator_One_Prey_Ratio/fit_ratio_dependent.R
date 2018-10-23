
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

general.NLL = function(attack, handling, exponent, initial, killed, predators, expttype, Pminus1, time=NULL){
	# DEBUG we should probably force some of these depending on model type
	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	attack <- exp(attack)
	handling <- exp(handling)
	exponent <- exp(exponent)

	if(is.null(time)){
		time <- 1
	}

	# expected number consumed
	Nconsumed <- general.pred(N0=initial, a=attack, h=handling, m=exponent, P=predators, T=time, expttype=expttype, Pminus1=Pminus1)

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

general.NLL.params = function(params, modeltype, initial, killed, predators, expttype, Pminus1, time=NULL){
	if(modeltype == "Holling I"){
		attack <- params[1]
		handling <- log(1E-7)
		exponent <- log(1E-7)
	}

	if(modeltype == "Holling II"){
		attack <- params[1]
		handling <- params[2]
		exponent <- log(1E-7)
	}
  
  if(modeltype == "Ratio I"){
    attack <- params[1]
    handling <- log(1E-7)
    exponent <- qlogis(1E-9)
  }
  
  if(modeltype == "Ratio II"){
		attack <- params[1]
		handling <- params[2]
		exponent <- qlogis(1E-9)
	}

	if(modeltype == "Hassell-Varley I"){
		attack <- params[1]
		handling <- log(1E-7)
		exponent <- params[2]

	}
	
	if(modeltype == "Hassell-Varley II"){
		attack <- params[1]
		handling <- params[2]
		exponent <- params[3]
	}

	nll <- general.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype, Pminus1=Pminus1)
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
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1)
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
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL)
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.ratioI.nloptr <- nloptr::nloptr(
  x0=x0.hl,
  eval_f=general.NLL.params,
  opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
  modeltype="Ratio I",
  initial=d$Nprey,
  killed=d$Nconsumed,
  predators=d$Npredator,
  time=d$Time,
  expttype=expttype,
  Pminus1=Pminus1
)

ffr.ratioI <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.ratioI.nloptr$solution[1]),
  fixed=list(
    handling = log(1E-7),
    exponent = log(1E-7)
  ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1)
  # optimizer="nlminb",
  # lower=c(attack=0, handling=0)
)

ffr.ratioII.nloptr <- nloptr::nloptr(
  x0=c(coef(ffr.ratioI)[1],log(1)),
  eval_f=general.NLL.params,
  opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
  modeltype="Ratio II",
  initial=d$Nprey,
  killed=d$Nconsumed,
  predators=d$Npredator,
  time=d$Time,
  expttype=expttype,
  Pminus1=Pminus1
)

ffr.ratioII <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.ratioII.nloptr$solution[1], handling=ffr.ratioII.nloptr$solution[2]),
  fixed=list(
    exponent = log(1E-7)
  ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL)
  # optimizer="nlminb",
  # lower=c(attack=0, handling=0)
)

ffr.HVI.nloptr <- nloptr::nloptr(
  x0=x0.hl,
  eval_f=general.NLL.params,
  opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
  modeltype="Hassel-Varley I",
  initial=d$Nprey,
  killed=d$Nconsumed,
  predators=d$Npredator,
  time=d$Time,
  expttype=expttype,
  Pminus1=Pminus1
)

ffr.HVI <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.HVI.nloptr$solution[1], exponent=ffr.HVI.nloptr$solution[2]),
  fixed=list(
    handling = log(1E-7)
  ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1)
  # optimizer="nlminb",
  # lower=c(attack=0, handling=0)
)

ffr.HVII.nloptr <- nloptr::nloptr(
  x0=c(coef(ffr.HVI)[1],log(1)),
  eval_f=general.NLL.params,
  opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
  modeltype="Hassel-Varley II",
  initial=d$Nprey,
  killed=d$Nconsumed,
  predators=d$Npredator,
  time=d$Time,
  expttype=expttype,
  Pminus1=Pminus1
)

ffr.HVII <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.HVII.nloptr$solution[1], handling=ffr.HVII.nloptr$solution[2], exponent=ffr.HVII.nloptr$solution[3]),
  # fixed=list(
  #   
  # ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL)
  # optimizer="nlminb",
  # lower=c(attack=0, handling=0)
)
