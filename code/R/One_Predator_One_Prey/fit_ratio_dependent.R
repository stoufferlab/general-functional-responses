
# libraries used below for the optimization/parameter estimation
library(bbmle)
library(nloptr)

# the likelihood functions are defined elsewhere stuff
source('../LogLikelihoods/ratio_dependent_single_predator_single_prey.R')

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