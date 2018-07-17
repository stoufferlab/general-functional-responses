
# EXTREME WEATHER WARNING:
# the fits produced by this version of the code can give better goodness-of-fit measures and lower AIC
# but could produce non-biologically-realistic estimates if extrapolated beyond the abundances/densities from the experiment

# libraries used below
library(bbmle)
library(nloptr)

# the likelihood functions are defined elsewhere
source('../LogLikelihoods/holling_like_one_predator_two_prey_unbounded.R')

####################################
registerDoParallel(cores=3)
####################################

# #####################################
# # series of functional response fits
# #####################################

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
		phi_ij = 0,
		phi_ji = 0
	),
	data=list(
		Ni = d$Prey.1.density,
		Nj = d$Prey.2.density,
		Ni_consumed = d$Prey.1.consumed,
		Nj_consumed = d$Prey.2.consumed,
		Npredators = d$Predators,
		expttype = "integrated"
	),
	skip.hessian = TRUE,
	# skip.hessian = FALSE,
	control=list(
		# trace=TRUE,
		# maxit=5000,
		# parscale=abs(coef(m2))
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
		phi_ij = 1,
		phi_ji = 1
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)

# fit a TYPE II functional response with phi_ji=1 and phi_ij a free parameter
ffr.typeII.ij.1 <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ij = 0.5
	),
	fixed=list(
		phi_ji = 1
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)

# fit a TYPE II functional response with phi_ji=0 and phi_ij a free parameter
ffr.typeII.ij.0 <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ij = 0.5
	),
	fixed=list(
		phi_ji = 0
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)

# fit a TYPE II functional response with phi_ij=1 and phi_ji a free parameter
ffr.typeII.1.ji <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ji = 0.5
	),
	fixed=list(
		phi_ij = 1
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)

# fit a TYPE II functional response with phi_ij=0 and phi_ji a free parameter
ffr.typeII.0.ji <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ji = 0.5
	),
	fixed=list(
		phi_ij = 0
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)

# fit a TYPE II functional response with phi_ij and phi_ji as free parameters
ffr.typeII.ij.ji <- bbmle::mle2(
	holling.like.1pred.2prey.NLL,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi_ji = 0.5,
		phi_ij = 0.5
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
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
		phi_ji = 0,
		phi_ij = 0
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)

# fit a TYPE II functional response with phi_ij = phi_ji and a single free parameter
ffr.typeII.ij.ij <- bbmle::mle2(
	holling.like.1pred.2prey.NLL.same_phi,
	start=list(
		attack_i = coef(ffr.typeI)["attack_i"],
		attack_j = coef(ffr.typeI)["attack_j"],
		handling_i = log(0.1),
		handling_j = log(0.1),
		phi = 0.5
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = TRUE
)


# fit a TYPE II functional response with phi_ij = phi_ji and a single free parameter
ffr.typeII.ij.ij <- bbmle::mle2(
	holling.like.1pred.2prey.NLL.same_phi,
	start=list(
		attack_i = coef(ffr.typeII.ij.ij)["attack_i"],
		attack_j = coef(ffr.typeII.ij.ij)["attack_j"],
		handling_i = coef(ffr.typeII.ij.ij)["handling_i"],
		handling_j = coef(ffr.typeII.ij.ij)["handling_j"],
		phi = coef(ffr.typeII.ij.ij)["phi"]
	),
	data=list(
		Ni=d$Prey.1.density,
		Nj=d$Prey.2.density,
		Ni_consumed=d$Prey.1.consumed,
		Nj_consumed=d$Prey.2.consumed,
		Npredators=d$Predators,
		expttype="integrated"
	),
	skip.hessian = FALSE,
	control=list(
		# trace=TRUE,
		# maxit=5000,
		parscale=abs(coef(ffr.typeII.ij.ij))
	)	
)



