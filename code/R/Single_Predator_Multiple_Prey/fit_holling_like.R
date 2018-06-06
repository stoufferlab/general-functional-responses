
# libraries used below
library(bbmle)
library(deSolve)
library(nloptr)

####################################
library(doParallel)
registerDoParallel(cores=3)
####################################

####################################
# generalizable functional responses
####################################

# predict the number of prey consumed in a replacement or integrated experiment
# for a single predator with multiple prey species and potential "handling-time exchangability" between those prey
general.pred = function(Ni, ai, hi, Nj, aj, hj, phi_ij, phi_ji, P, T, expttype=c("integrated","replacement")){
	expttype <- match.arg(expttype)

	print(c(ai=ai, aj=aj, hi=hi, hj=hj, phi_ij=phi_ij, phi_ji=phi_ji))

	# continuous replacement of consumed prey
	if(expttype=="replacement"){
		# there is a common demoninator for both species
		denom <- (1 + ai * hi * Ni) * (1 + aj * hj * Nj) - phi_ij * phi_ji * ai * hi * Ni * aj * hj * Nj

		# predicted numerator and number consumed for species i
		numeri <- ai * Ni * (1 + (1 - phi_ij) * aj * hj * Nj)
		Ni_e <- (numeri / denom) * P * T

		# predicted numerator and number consumed of species j
		numerj <- aj * Nj * (1 + (1 - phi_ji) * ai * hi * Ni)
		Nj_e <- (numerj / denom) * P * T

		# combined predicted number of prey consumed
		Ne <- cbind(Ni_e, Nj_e)
	}

	# non-replacement of consumed prey
	if(expttype=="integrated"){
		# the derivative in the integration
		func <- function(t, y, parms){
			Ni <- y[1]
			Nj <- y[2]
			with(as.list(parms),
			{
				denom <- (1 + ai * hi * Ni) * (1 + aj * hj * Nj) - phi_ij * phi_ji * ai * hi * Ni * aj * hj * Nj

				numeri <- ai * Ni * (1 + (1 - phi_ij) * aj * hj * Nj)
				numerj <- aj * Nj * (1 + (1 - phi_ji) * ai * hi * Ni)

				dNidt <- -1 * (numeri / denom) * P
				dNjdt <- -1 * (numerj / denom) * P

				return(list(c(dNidt, dNjdt)))
			})
		}
		
		# container for the model output
		# Ne <- matrix(NA, length(Ni), 2)

		# DEBUG this can be broken down into replicates since the prediction of the model is the same for the same 'treatment' conditions
		Ne <- foreach(i=1:length(Ni), .combine=rbind) %dopar% {
			# initial conditions
			N0 <- c(Ni[i], Nj[i])

			# list all parameters like ode expects them
			parms <- c(ai=ai, aj=aj, hi=hi, hj=hj, phi_ij=phi_ij, phi_ji=phi_ji, P=P[i])
		
			# integrate 1000 time steps from 0 to the duration of experiment
			times <- seq(0,T[i],length.out=100)
		
			# integrate and determine the solution
			Nf <- deSolve::ode(N0, times, func, parms)
			Nf <- Nf[nrow(Nf),2:3]

			# the number consumed is the difference between what we started with and what is left
			# Ne[i,] <- N0 - Nf
			N0 - Nf
		}
		# print(parms[1:6])
	}

	return(Ne)
}

general.NLL = function(attack_i, handling_i, attack_j, handling_j, phi_ij, phi_ji, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time=rep(1,length(Ni))){
	# DEBUG we should probably force some of these here depending on model type
	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	# standard Holling Type II parameters
	attack_i <- exp(attack_i)
	attack_j <- exp(attack_j)
	handling_i <- exp(handling_i)
	handling_j <- exp(handling_j)

	# new parameters	
	phi_ij <- plogis(phi_ij)
	phi_ji <- plogis(phi_ji)

	# expected number consumed
	Nconsumed <- general.pred(
		Ni=Ni,
		Nj=Nj,
		ai=attack_i,
		aj=attack_j,
		hi=handling_i,
		hj=handling_j,
		phi_ij=phi_ij,
		phi_ji=phi_ji,
		P=Npredators,
		T=time,
		expttype=expttype
	)

	# negative log likelihood based on proportion consumed (no replacement)
	if(expttype=="integrated"){
		nll <- -sum(dbinom(Ni_consumed, prob=pmax(0,Nconsumed[,1]/Ni,na.rm=TRUE), size=Ni, log=TRUE))
		nll <- nll - sum(dbinom(Nj_consumed, prob=pmax(0,Nconsumed[,2]/Nj,na.rm=TRUE), size=Nj, log=TRUE))
	}

	# negative log likelihood based on total number consumed (replacement)
	if(expttype=="replacement"){
		nll <- -sum(dpois(Ni_consumed, Nconsumed[,1], log=TRUE))
		nll <- nll - sum(dpois(Nj_consumed, Nconsumed[,2], log=TRUE))
	}

	return(nll)
}

general.NLL.same_phi = function(attack_i, handling_i, attack_j, handling_j, phi, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time=rep(1,length(Ni))){
	nll <- general.NLL(attack_i, handling_i, attack_j, handling_j, phi, phi, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time)
	return(nll)
}

# general.NLL.params = function(params, modeltype, initial, killed, predators, time, expttype, Pminus1){
# 	if(modeltype == "Holling I"){
# 		attack <- params[1]
# 		handling <- -20
# 		interference <- -20
# 		phi_numer <- 20
# 		phi_denom <- 20	
# 	}

# 	if(modeltype == "Holling II"){
# 		attack <- params[1]
# 		handling <- params[2]
# 		interference <- -20
# 		phi_numer <- 20
# 		phi_denom <- 20	
# 	}

# 	if(modeltype == "Beddington-DeAngelis"){
# 		attack <- params[1]
# 		handling <- params[2]
# 		interference <- params[3]
# 		phi_numer <- 20
# 		phi_denom <- 20	
# 	}

# 	if(modeltype == "Crowley-Martin"){
# 		attack <- params[1]
# 		handling <- params[2]
# 		interference <- params[3]
# 		phi_numer <- 20
# 		phi_denom <- -20
# 	}

# 	if(modeltype == "Stouffer-Novak I"){
# 		attack <- params[1]
# 		handling <- params[2]
# 		interference <- params[3]
# 		phi_numer <- 20
# 		phi_denom <- params[4]
# 	}

# 	if(modeltype == "Stouffer-Novak II"){
# 		attack <- params[1]
# 		handling <- params[2]
# 		interference <- params[3]
# 		phi_numer <- params[4]
# 		phi_denom <- 20
# 	}
	
# 	if(modeltype == "Stouffer-Novak III"){
# 		attack <- params[1]
# 		handling <- params[2]
# 		interference <- params[3]
# 		phi_numer <- params[4]
# 		phi_denom <- params[5]
# 	}

# 	nll <- general.NLL(attack=attack, handling=handling, interference=interference, phi_numer=phi_numer, phi_denom=phi_denom, initial=initial, killed=killed, predators=predators, time=time, expttype=expttype, Pminus1=Pminus1)
# 	return(nll)
# }

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
d <- read.csv("test_data.csv")
d <- d[rowSums(is.na(d))==0,]

# fit a TYPE I functional response treating the predator as specialized on both prey
ffr.typeI <- bbmle::mle2(
	general.NLL,
	start=list(
		attack_i = log(0.1),
		attack_j = log(0.1)
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
ffr.typeII.specialist <- bbmle::mle2(
	general.NLL,
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

# fit a TYPE II functional response treating the predator as specialized on both prey
ffr.typeII.ij <- bbmle::mle2(
	general.NLL,
	start=list(
		attack_i = coef(ffr.typeII.specialist)["attack_i"],
		attack_j = coef(ffr.typeII.specialist)["attack_j"],
		handling_i = coef(ffr.typeII.specialist)["handling_i"],
		handling_j = coef(ffr.typeII.specialist)["handling_j"],
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

# fit a TYPE II functional response treating the predator as specialized on both prey
ffr.typeII.ji <- bbmle::mle2(
	general.NLL,
	start=list(
		attack_i = coef(ffr.typeII.specialist)["attack_i"],
		attack_j = coef(ffr.typeII.specialist)["attack_j"],
		handling_i = coef(ffr.typeII.specialist)["handling_i"],
		handling_j = coef(ffr.typeII.specialist)["handling_j"],
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

# fit a TYPE II functional response treating the predator as specialized on both prey
ffr.typeII.both <- bbmle::mle2(
	general.NLL,
	start=list(
		attack_i = coef(ffr.typeII.specialist)["attack_i"],
		attack_j = coef(ffr.typeII.specialist)["attack_j"],
		handling_i = coef(ffr.typeII.specialist)["handling_i"],
		handling_j = coef(ffr.typeII.specialist)["handling_j"],
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

# fit a TYPE II functional response treating the predator as specialized on both prey
ffr.typeII.generalist <- bbmle::mle2(
	general.NLL,
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

# fit a TYPE II functional response treating the predator as specialized on both prey
ffr.typeII.same_phi <- bbmle::mle2(
	general.NLL.same_phi,
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

