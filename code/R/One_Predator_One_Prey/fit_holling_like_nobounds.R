
# EXTREME WEATHER WARNING:
# the fits produced by this version of the code can give better goodness-of-fit measures and lower AIC
# but could produce non-biologically-realistic estimates if extrapolated beyond the abundances/densities from the experiment

# libraries used below
library(bbmle)
library(nloptr)

# the likelihood functions are defined elsewhere
source('../LogLikelihoods/holling_like_one_predator_one_prey_unbounded.R')

#####################################
# series of functional response fits
#####################################

ffr.hollingI.nloptr <- nloptr::nloptr(
	x0=x0.hl,
	eval_f=holling.like.1pred.1prey.NLL,
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
	holling.like.1pred.1prey.NLL,
	vecpar=TRUE,
	start=list(attack=ffr.hollingI.nloptr$solution[1]),
	fixed=list(
		handling = 0,
		interference = 0,
		phi_numer = 1,
		phi_denom = 1
	),
	data=list(
		modeltype="Holling I",
		initial=d$Nprey,
		killed=d$Nconsumed,
		predators=d$Npredator,
		time=d$Time,
		expttype=expttype,
		Pminus1=Pminus1
	)
	#skip.hessian = TRUE
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

ffr.hollingII.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingI)[1],log(1)),
	eval_f=holling.like.1pred.1prey.NLL,
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
	holling.like.1pred.1prey.NLL,
	vecpar=TRUE,
	start=list(
		attack=ffr.hollingII.nloptr$solution[1],
		handling=ffr.hollingII.nloptr$solution[2]
	),
	fixed=list(
		interference = 0,
		phi_numer = 1,
		phi_denom = 1
	),
	data=list(modeltype="Holling II", initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1)
	#skip.hessian = TRUE
	# optimizer="nlminb",
	# lower=c(attack=0, handling=0)
)

# ffr.bd.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],log(1)),
# 	eval_f=holling.like.1pred.1prey.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	modeltype="Beddington-DeAngelis",
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# ffr.bd <- bbmle::mle2(
# 	holling.like.1pred.1prey.NLL,
# 	start=list(attack=ffr.bd.nloptr$solution[1], handling=ffr.bd.nloptr$solution[2], interference=ffr.bd.nloptr$solution[3]),
# 	fixed=list(
# 		phi_numer = 1,
# 		phi_denom = 1
# 	),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
# 	#skip.hessian = TRUE
# )

# ffr.cm.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],log(1)),
# 	eval_f=holling.like.1pred.1prey.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	modeltype="Crowley-Martin",
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# ffr.cm <- bbmle::mle2(
# 	holling.like.1pred.1prey.NLL,
# 	start=list(attack=ffr.cm.nloptr$solution[1], handling=ffr.cm.nloptr$solution[2], interference=ffr.cm.nloptr$solution[3]),
# 	fixed=list(
# 		phi_numer = 1,
# 		phi_denom = 0
# 	),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL),
# 	#skip.hessian = TRUE
# 	# optimizer="nlminb",
# 	# lower=c(attack=0, handling=0)
# )

# optimize with nloptr

ffr.sn1.nloptr <- nloptr::nloptr(
	x0=c(coef(ffr.hollingII)[1:2],log(1),0.5),
	eval_f=holling.like.1pred.1prey.NLL,
	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
	modeltype="Stouffer-Novak I",
	initial=d$Nprey,
	killed=d$Nconsumed,
	predators=d$Npredator,
	time=d$Time,
	expttype=expttype,
	Pminus1=Pminus1
)

# WARNING THIS IS FOR TESTING ONLY
# ffr.sn1.nloptr$solution <- c(coef(ffr.hollingII)[1:2],log(1),0.5)

# use mle2 because it has handy machinery

ffr.sn1 <- bbmle::mle2(
	holling.like.1pred.1prey.NLL,
	vecpar=TRUE,
	start=list(
		attack=ffr.sn1.nloptr$solution[1],
		handling=ffr.sn1.nloptr$solution[2],
		interference=ffr.sn1.nloptr$solution[3],
		phi_denom=ffr.sn1.nloptr$solution[4]
	),
	fixed=list(
		phi_numer = 1
	),
	data=list(
		modeltype="Stouffer-Novak I",
		initial=d$Nprey,
		killed=d$Nconsumed,
		predators=d$Npredator,
		time=d$Time,
		expttype=expttype,
		Pminus1=Pminus1
	),
	skip.hessian = TRUE
)

# try to get closer with mle2 by using the parscale argument, and keep the hessian for confidence intervals

ffr.sn1 <- bbmle::mle2(
	holling.like.1pred.1prey.NLL,
	vecpar=TRUE,
	start=list(
		attack=coef(ffr.sn1)[1],
		handling=coef(ffr.sn1)[2],
		interference=coef(ffr.sn1)[3],
		phi_denom=coef(ffr.sn1)[5]
	),
	fixed=list(
		phi_numer = 1
	),
	data=list(
		modeltype="Stouffer-Novak I",
		initial=d$Nprey,
		killed=d$Nconsumed,
		predators=d$Npredator,
		time=d$Time,
		expttype=expttype,
		Pminus1=Pminus1
	),
	control=list(
		# trace=TRUE,
		# maxit=5000,
		parscale=abs(coef(ffr.sn1)[c(1:3,5)])
	)
)


# # DEBUG keeps failing for some datasets
# ffr.sn2.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],log(1),0.5),
# 	eval_f=holling.like.1pred.1prey.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# ffr.sn2 <- bbmle::mle2(
# 	holling.like.1pred.1prey.NLL,
# 	start=list(attack=ffr.sn2.nloptr$solution[1], handling=ffr.sn2.nloptr$solution[2], interference=ffr.sn2.nloptr$solution[3], phi_numer=ffr.sn2.nloptr$solution[4]),
# 	fixed=list(
# 		phi_denom=1
# 	),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
# 	#skip.hessian = TRUE
# )

# ffr.sn3.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],log(1),0.5,0.5),
# 	eval_f=holling.like.1pred.1prey.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	modeltype="Stouffer-Novak III",
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# ffr.sn3 <- bbmle::mle2(
# 	holling.like.1pred.1prey.NLL,
# 	start=list(attack=ffr.sn3.nloptr$solution[1], handling=ffr.sn3.nloptr$solution[2], interference=ffr.sn3.nloptr$solution[3], phi_numer=ffr.sn3.nloptr$solution[4], phi_denom=ffr.sn3.nloptr$solution[5]),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
# 	#skip.hessian = TRUE
# )
