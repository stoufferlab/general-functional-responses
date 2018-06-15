
# EXTREME WEATHER WARNING:
# the fits produced by this version of the code can give better goodness-of-fit measures and lower AIC
# but could produce non-biologically-realistic estimates if extrapolated beyond the abundances/densities from the experiment

# libraries used below
library(bbmle)
library(nloptr)

# the likelihood functions are defined elsewhere
source('../LogLikelihoods/holling_like_one_predator_two_prey.R')

####################################
registerDoParallel(cores=3)
####################################

# #####################################
# # series of functional response fits
# #####################################

# ffr.hollingI.nloptr <- nloptr::nloptr(
# 	x0=x0.hl,
# 	eval_f=general.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	modeltype="Holling I",
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# read in some test data for the time being
d <- read.csv("test_data_1.csv")
d <- d[rowSums(is.na(d))==0,]
# d <- read.csv("test_data_2.csv")
# d$Prey.1.consumed <- round(d$Prey.1.consumed+rnorm(nrow(d),0,0.001))

# fit a TYPE I functional response treating the predator as specialized on both prey
ffr.typeI <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = log(0.2),
		attack_j = log(0.2)
	),
	fixed=list(
		handling_i = log(1E-7),
		handling_j = log(1E-7),
		phi_ij = qlogis(1E-9),
		phi_ji = qlogis(1E-9)),
	data=list(
		Ni = d$Prey.1.density,
		Nj = d$Prey.2.density,
		Ni_consumed = d$Prey.1.consumed,
		Nj_consumed = d$Prey.2.consumed,
		Npredators = d$Predators,
		expttype = "integrated"
	)
)

# fit a TYPE II functional response treating the predator as specialized on both prey
ffr.typeII.1.1 <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1)
	),
	fixed=list(
		phi_ij=qlogis(1E-9),
		phi_ji=qlogis(1E-9)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ji=1 and phi_ij a free parameter
ffr.typeII.ij.1 <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ij = qlogis(0.5)
	),
	fixed=list(
		phi_ji=qlogis(1E-9)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ji=0 and phi_ij a free parameter
ffr.typeII.ij.0 <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ij = qlogis(0.5)
	),
	fixed=list(
		phi_ji=-1 * qlogis(1E-9)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ij=1 and phi_ji a free parameter
ffr.typeII.1.ji <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ji = qlogis(0.5)
	),
	fixed=list(
		phi_ij=qlogis(1E-9)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ij=0 and phi_ji a free parameter
ffr.typeII.0.ji <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ji = qlogis(0.5)
	),
	fixed=list(
		phi_ij = -1 * qlogis(1E-9)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ij and phi_ji as free parameters
ffr.typeII.ij.ji <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ji = qlogis(0.5),
		phi_ij=qlogis(0.5)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ij=0 and phi_ji=0 as free parameters
ffr.typeII.0.0 <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1)
	),
	fixed=list(
		phi_ji = -1 * qlogis(1E-9),
		phi_ij = -1 * qlogis(1E-9)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)

# fit a TYPE II functional response with phi_ij = phi_ji and a single free parameter
ffr.typeII.ij.ij <- bbmle::mle2(
	holling.like.1pred.2prey.NLL.same_phi,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi = qlogis(0.5)
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	)
)


# ffr.hollingII.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingI)[1],0),
# 	eval_f=general.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	modeltype="Holling II",
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# ffr.hollingII <- bbmle::mle2(
# 	general.NLL,
# 	start=list(attack=ffr.hollingII.nloptr$solution[1], handling=ffr.hollingII.nloptr$solution[2]),
# 	fixed=list(interference=-20, phi_numer=20, phi_denom=20),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL),
# 	# optimizer="nlminb",
# 	# lower=c(attack=0, handling=0)
# )

# ffr.bd.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],0),
# 	eval_f=general.NLL.params,
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
# 	general.NLL,
# 	start=list(attack=ffr.bd.nloptr$solution[1], handling=ffr.bd.nloptr$solution[2], interference=ffr.bd.nloptr$solution[3]),
# 	fixed=list(phi_numer=20, phi_denom=20),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
# )

# ffr.cm.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],0),
# 	eval_f=general.NLL.params,
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
# 	general.NLL,
# 	start=list(attack=ffr.cm.nloptr$solution[1], handling=ffr.cm.nloptr$solution[2], interference=ffr.cm.nloptr$solution[3]),
# 	fixed=list(phi_numer=20, phi_denom=-20),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1, params=NULL),
# 	# optimizer="nlminb",
# 	# lower=c(attack=0, handling=0)
# )

# ffr.sn1.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],0,0),
# 	eval_f=general.NLL.params,
# 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# 	modeltype="Stouffer-Novak I",
# 	initial=d$Nprey,
# 	killed=d$Nconsumed,
# 	predators=d$Npredator,
# 	time=d$Time,
# 	expttype=expttype,
# 	Pminus1=Pminus1
# )

# ffr.sn1 <- bbmle::mle2(
# 	general.NLL,
# 	start=list(attack=ffr.sn1.nloptr$solution[1], handling=ffr.sn1.nloptr$solution[2], interference=ffr.sn1.nloptr$solution[3], phi_denom=ffr.sn1.nloptr$solution[4]),
# 	fixed=list(phi_numer=20),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
# )

# # DEBUG keeps failing for some datasets
# # ffr.sn2.nloptr <- nloptr::nloptr(
# # 	x0=c(coef(ffr.hollingII)[1:2],0,0),
# # 	eval_f=sn1.NLL,
# # 	opts=list(print_level=0, algorithm='NLOPT_LN_SBPLX', maxeval=1E5),
# # 	initial=d$Nprey,
# # 	killed=d$Nconsumed,
# # 	predators=d$Npredator,
# # 	time=d$Time,
# # 	expttype=expttype,
# # )

# # ffr.sn2 <- bbmle::mle2(
# # 	general.NLL,
# # 	start=list(attack=ffr.sn2.nloptr$solution[1], handling=ffr.sn2.nloptr$solution[2], interference=ffr.sn2.nloptr$solution[3], phi_numer=ffr.sn2.nloptr$solution[4]),
# # 	fixed=list(phi_denom=20),
# # 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype),
# # )

# ffr.sn3.nloptr <- nloptr::nloptr(
# 	x0=c(coef(ffr.hollingII)[1:2],0,0,0),
# 	eval_f=general.NLL.params,
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
# 	general.NLL,
# 	start=list(attack=ffr.sn3.nloptr$solution[1], handling=ffr.sn3.nloptr$solution[2], interference=ffr.sn3.nloptr$solution[3], phi_numer=ffr.sn3.nloptr$solution[4], phi_denom=ffr.sn3.nloptr$solution[5]),
# 	data=list(initial=d$Nprey, killed=d$Nconsumed, predators=d$Npredator, time=d$Time, expttype=expttype, Pminus1=Pminus1),
# )

