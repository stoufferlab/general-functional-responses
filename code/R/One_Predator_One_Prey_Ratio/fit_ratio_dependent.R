
# libraries used below
library(bbmle)
library(lamW)
library(nloptr)

###############################
# generalizable functional responses
###############################

general.pred = function(N0, a, h, m, P, T, expttype=c("integrated","replacement")){
	expttype <- match.arg(expttype)

	# In contrast to Holling-type, P_interfering is always P for ratio models
	
	if(expttype=="integrated"){
	  # h <- h * (1 + (1 - phi_numer * phi_denom) * a * h * c * P_interfering)
	  # Q <- (1 + c * P_interfering)
	  # X <- (1 + (1 - phi_numer) * c * P_interfering)
	  # N <- N0 - (Q / (a * h)) * lamW::lambertW0(((a * h * N0)/ Q) * exp(- (a / Q) * (X * P * T - h * N0)))
	}

	if(expttype=="replacement"){
		# numer <- a * N0 * (P ^ -m)
		# denom <- 1 + a * h * N0 * (P ^ -m)
		# N <- (numer / denom) * P * T
	}
	
	return(N)
}

general.NLL = function(attack, handling, exponent, initial, killed, predators, expttype, time=NULL){
  
	attack <- exp(attack)
	handling <- exp(handling)
	exponent <- exp(exponent)

	if(is.null(time)){
		time <- 1
	}

	# expected number consumed
	Nconsumed <- general.pred(N0=initial, a=attack, h=handling, m=exponent, P=predators, T=time, expttype=expttype)

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

general.NLL.params = function(params, modeltype, initial, killed, predators, expttype, time=NULL){
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
    exponent <- 0
  }
  
  if(modeltype == "Ratio II"){
		attack <- params[1]
		handling <- params[2]
		exponent <- 0
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

	nll <- general.NLL(attack=attack, handling=handling, exponent=exponent, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype)
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
	expttype=expttype
)

ffr.hollingI <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.hollingI.nloptr$solution[1]),
	fixed=list(
		handling = log(1E-7),
		exponent = log(1E-7)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype)
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
	expttype=expttype
)

ffr.hollingII <- bbmle::mle2(
	general.NLL,
	start=list(attack=ffr.hollingII.nloptr$solution[1], handling=ffr.hollingII.nloptr$solution[2]),
	fixed=list(
		exponent = log(1E-7)
	),
	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, params=NULL)
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
  expttype=expttype
)

ffr.ratioI <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.ratioI.nloptr$solution[1]),
  fixed=list(
    handling = log(1E-7),
    exponent = 0
  ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype)
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
  expttype=expttype
)

ffr.ratioII <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.ratioII.nloptr$solution[1], handling=ffr.ratioII.nloptr$solution[2]),
  fixed=list(
    exponent = 0
  ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, params=NULL)
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
  expttype=expttype
)

ffr.HVI <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.HVI.nloptr$solution[1], exponent=ffr.HVI.nloptr$solution[2]),
  fixed=list(
    handling = log(1E-7)
  ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype)
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
  expttype=expttype
)

ffr.HVII <- bbmle::mle2(
  general.NLL,
  start=list(attack=ffr.HVII.nloptr$solution[1], handling=ffr.HVII.nloptr$solution[2], exponent=ffr.HVII.nloptr$solution[3]),
  # fixed=list(
  #   
  # ),
  data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, params=NULL)
)
