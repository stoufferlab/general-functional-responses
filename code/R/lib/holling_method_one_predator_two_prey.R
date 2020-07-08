
#############################################
# holling-like functional responses
#############################################

# libraries used below
library(bbmle)
library(nloptr)
library(odeintr)
# library(doParallel); registerDoParallel(cores=1)

# define the ode in C++ format to use odeintr
holling.like.1pred.2prey.sys = '
	// hybrid functional response for one predator two prey
	dxdt[0] = -P * (ai * x[0] * (1 + (1 - phi_ij) * aj * hj * x[1])) / ((1 + ai * hi * x[0]) * (1 + aj * hj * x[1]) - phi_ij * phi_ji * ai * hi * x[0] * aj * hj * x[1]);
	
	// consumption rate cannot be positive
	if(dxdt[0] > 0) dxdt[0] = 0;

	// hybrid functional response for one predator two prey
	dxdt[1] = -P * (aj * x[1] * (1 + (1 - phi_ji) * ai * hi * x[0])) / ((1 + ai * hi * x[0]) * (1 + aj * hj * x[1]) - phi_ij * phi_ji * ai * hi * x[0] * aj * hj * x[1]);
	
	// consumption rate cannot be positive
	if(dxdt[1] > 0) dxdt[1] = 0;
'

# compile the above C++ code into something we can run in R
odeintr::compile_sys(
	"hl_1pred_2prey",
	holling.like.1pred.2prey.sys,
	pars = c("ai", "aj", "hi", "hj", "phi_ij", "phi_ji", "P") #,
	# method = "bsd"
)

# predict the number of prey consumed in a replacement or integrated experiment
# for a single predator with multiple prey species and potential "handling-time exchangability" between those prey
holling.like.1pred.2prey = function(ai, hi, aj, hj, phi_ij, phi_ji, Ni, Nj, P, T, replacement){
	# print(c(ai=ai, aj=aj, hi=hi, hj=hj, phi_ij=phi_ij, phi_ji=phi_ji))

	# continuous replacement of consumed prey
	if(replacement){
		# there is a common demoninator for both species
		denom <- (1 + ai * hi * Ni) * (1 + aj * hj * Nj) - phi_ij * phi_ji * ai * hi * Ni * aj * hj * Nj

		# predicted numerator and number consumed for species i
		numeri <- ai * Ni * (1 + (1 - phi_ij) * aj * hj * Nj)
		
		# predicted number consumed of species i
		Ni_e <- (numeri / denom) * P * T

		# make sure we do not predict negative consumption 
		Ni_e <- pmax(0, Ni_e)

		# predicted numerator and number consumed of species j
		numerj <- aj * Nj * (1 + (1 - phi_ji) * ai * hi * Ni)
		
		# predicted number consumed of species j
		Nj_e <- (numerj / denom) * P * T

		# make sure we do not predict negative consumption
		Nj_e <- pmax(0, Nj_e)

		# combined predicted number of prey consumed across both prey species
		Neaten <- cbind(Ni_e, Nj_e)

		return(Neaten)
	}

	# non-replacement of consumed prey
	if(!replacement){
		# DEBUG to speed up this can be broken down into identical replicates based on initial abundances since the prediction of the model is the same for the same 'treatment' conditions
		Neaten <- matrix(NA, nrow=length(Ni), ncol=2)
		for(i in seq.int(length(Ni))){
			# initial conditions
			N0 <- c(Ni[i], Nj[i])

			# set parameters within ode solver
			hl_1pred_2prey_set_params(ai=ai, aj=aj, hi=hi, hj=hj, phi_ij=phi_ij, phi_ji=phi_ji, P=P[i])

			# calculate the final number of prey integrating the odes
			# note that the dynamic model (defined above) does not permit "negative" consumption
			Nfinal <- hl_1pred_2prey(N0, T[i], T[i]/1000.)

			# we technically only need the last row since this is the final "abundance"
			Nfinal <- as.numeric(Nfinal[nrow(Nfinal),2:3])

			# on rare occasions we dip into negative final abundances (that tend to be extremely small)
			Nfinal <- pmax(0, Nfinal)

			# let's also make sure no prey magically appear
			Nfinal <- pmin(Nfinal, N0)

			# the number consumed is the difference between what we started with and what is left
			Neaten[i,] <- N0 - Nfinal
		}

		return(Neaten)
	}

	stop()
}

holling.like.1pred.2prey.predict = function(params, Ni, Nj, Npredators, replacement, modeltype, phi.transform, time){
	if(modeltype=="Holling I"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- 0
		handling_j <- 0
		phi_ij <- 0
		phi_ji <- 0
	}

	if(modeltype=="Holling II Specialist Specialist"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- 0
		phi_ji <- 0
	}

	if(modeltype=="Holling II Specialist Generalist"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- 0
		phi_ji <- 1
	}

	if(modeltype=="Holling II Generalist Specialist"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- 1
		phi_ji <- 0
	}

	if(modeltype=="Holling II Generalist Generalist"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- 1
		phi_ji <- 1
	}

	if(modeltype=="Holling II Specialist Hybrid"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- 0
		phi_ji <- phi.transform(params[5])
	}

	if(modeltype=="Holling II Generalist Hybrid"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- 1
		phi_ji <- phi.transform(params[5])
	}		

	if(modeltype=="Holling II Hybrid Specialist"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- phi.transform(params[5])
		phi_ji <- 0
	}

	if(modeltype=="Holling II Hybrid Generalist"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- phi.transform(params[5])
		phi_ji <- 1
	}
	
	if(modeltype=="Holling II Hybrid Hybrid"){
		attack_i <- exp(params[1])
		attack_j <- exp(params[2])
		handling_i <- exp(params[3])
		handling_j <- exp(params[4])
		phi_ij <- phi.transform(params[5])
		phi_ji <- phi.transform(params[6])
	}

	# if no times are specified then normalize to time=1
	if(is.null(time)){
		time <- rep(1,length(Ni))
	}

	# expected number consumed
	Nconsumed <- holling.like.1pred.2prey(
		ai=attack_i,
		aj=attack_j,
		hi=handling_i,
		hj=handling_j,
		phi_ij=phi_ij,
		phi_ji=phi_ji,
		Ni=Ni,
		Nj=Nj,
		P=Npredators,
		T=time,
		replacement=replacement
	)

	return(Nconsumed)
}

holling.like.1pred.2prey.NLL = function(params, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, replacement, modeltype, phi.transform, time=NULL){
	# # speed things up by identifying equivalent treatments
	# d <- data.frame(Ni=Ni, Nj=Nj, Npredators=Npredators)

	# # find unique rows and extract their treatment levels
	# d.uniq <- unique(d)

	# # get the predictions
	# Nconsumed <- holling.like.1pred.2prey.predict(
	# 	params,
	# 	d.uniq$Ni,
	# 	d.uniq$Nj,
	# 	d.uniq$Npredators,
	# 	replacement,
	# 	modeltype,
	# 	phi.transform,
	# 	time
	# )

	# # re-expand out to the full dataset
	# d.uniq$Ni_consumed <- Nconsumed[,1]
	# d.uniq$Nj_consumed <- Nconsumed[,2]

	# print(d.uniq)

	# d <- merge(d, d.uniq, by=c("Ni","Nj","Npredators"))

	# print(d)

	# # re-extract "full" version as Nconsumed
	# Nconsumed <- d[,c("Ni_consumed","Nj_consumed")]

	# print(cbind(Nconsumed,Ni_consumed,Nj_consumed))

	# get the predictions
	Nconsumed <- holling.like.1pred.2prey.predict(
		params,
		Ni,
		Nj,
		Npredators,
		replacement,
		modeltype,
		phi.transform,
		time
	)

	# print(dim(Nconsumed))

	# XXX

	# DEBUG check whether the parameters give biological valid predictions
	# should this occur here or above?

	# disallow biologically implausible predictions
	if(any(Nconsumed < 0) || any(!is.finite(as.matrix(Nconsumed)))){
		return(Inf)
	}

	# negative log likelihood based on total number consumed (replacement)
	if(replacement){
		nll <- -sum(dpois(Ni_consumed, Nconsumed[,1], log=TRUE))
		nll <- nll - sum(dpois(Nj_consumed, Nconsumed[,2], log=TRUE))
		return(nll)
	}

	# negative log likelihood based on proportion consumed (no replacement)
	if(!replacement){
		nll <- -sum(dbinom(Ni_consumed, prob=pmax(0, ifelse(Ni == 0, 0, Nconsumed[,1] / Ni)), size=Ni, log=TRUE))
		nll <- nll - sum(dbinom(Nj_consumed, prob=pmax(0, ifelse(Nj == 0, 0, Nconsumed[,2] / Nj)), size=Nj, log=TRUE))
		return(nll)
	}

	# if we made it this far something is wrong
	stop()
}

# holling.like.1pred.2prey.NLL.same_phi = function(attack_i, handling_i, attack_j, handling_j, phi, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time=rep(1,length(Ni))){
# 	nll <- holling.like.1pred.2prey.NLL(attack_i, handling_i, attack_j, handling_j, phi, phi, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time)
# 	return(nll)
# }

parnames(holling.like.1pred.2prey.NLL) <- c(
	'attack_i',
	'attack_j',
	'handling_i',
	'handling_j',
	'phi_ij',
	'phi_ji'
)

# DEBUG: should allow starting parameters to be specified since fitting is very slow for this approach
# given data (d), study info (s), and modeltype (e.g., "Holling I"), fit functional response data
fit.holling.like <- function(d, s, modeltype, phi.transform=identity, nloptr.control=list(), mle2.control=list(), ...){
	# estimate starting values from the data using linear regression
	start <- c(
		log(coef(lm(d$Nconsumed1~0+I(d$Npredator * d$Nprey1)))),
		log(coef(lm(d$Nconsumed2~0+I(d$Npredator * d$Nprey2))))
	)
	names(start) <- c("attack_i", "attack_j")

	# fit Holling Type I via MLE with above starting parameter value
	hollingI.via.sbplx <- nloptr::sbplx(
		x0 = start,
		fn = holling.like.1pred.2prey.NLL,
		Ni = d$Nprey1,
		Nj = d$Nprey2,
		Ni_consumed = d$Nconsumed1,
		Nj_consumed = d$Nconsumed2,
		Npredators = d$Npredator,
		time=d$Time,
		replacement = s$replacement,
		modeltype = "Holling I",
		control = nloptr.control #,
		# ...
	)

	# print(hollingI.via.sbplx$value)

	mle2.start <- as.list(hollingI.via.sbplx$par)
	names(mle2.start) <- names(start)

	# fit a purely linear model to get appropriate starting values
	hollingI.via.mle2 <- bbmle::mle2(
		holling.like.1pred.2prey.NLL,
		start = mle2.start,
		data = list(
			Ni = d$Nprey1,
			Nj = d$Nprey2,
			Ni_consumed = d$Nconsumed1,
			Nj_consumed = d$Nconsumed2,
			Npredators = d$Npredator,
			time=d$Time,
			replacement = s$replacement,
			modeltype = "Holling I"
		),
		vecpar = TRUE,
		control = mle2.control,
		# ...
	)

	if(modeltype == "Holling I"){
		# print(hollingI.via.mle2@coef)
		return(hollingI.via.mle2)
	}else{
		# if moving to a more complex model, fit a specialist-specialist holling type II as the next starting point
		start <- list(
			attack_i = coef(hollingI.via.mle2)["attack_i"],
			attack_j = coef(hollingI.via.mle2)["attack_j"],
			handling_i = log(1),
			handling_j = log(1)
		)

		# DEBUG sometimes Holling Type I doesn't provide the best starting guess for the more complex models
		# DEBUG we should think about whether or not there are better options to make sure we can be confident about the final fits

		# # fit Holling Type II Generalist Generalist with above starting parameter values
		# hollingII.via.sbplx <- nloptr::sbplx(
		# 	x0 = unlist(start),
		# 	fn = holling.like.1pred.2prey.NLL,
		# 	Ni = d$Nprey1,
		# 	Nj = d$Nprey2,
		# 	Ni_consumed = d$Nconsumed1,
		# 	Nj_consumed = d$Nconsumed2,
		# 	Npredators = d$Npredator,
		# 	replacement = s$replacement,
		# 	modeltype = "Holling II Generalist Generalist",
		# 	control = nloptr.control #,
		# 	# ...
		# )

		# mle2.start <- as.list(hollingII.via.sbplx$par)
		# names(mle2.start) <- names(start)

		# # refit this model with mle2
		# hollingII.via.mle2 <- bbmle::mle2(
		# 	holling.like.1pred.2prey.NLL,
		# 	start=mle2.start,
		# 	data=list(
		# 		Ni = d$Nprey1,
		# 		Nj = d$Nprey2,
		# 		Ni_consumed = d$Nconsumed1,
		# 		Nj_consumed = d$Nconsumed2,
		# 		Npredators = d$Npredator,
		# 		replacement = s$replacement,
		# 		modeltype = "Holling II Generalist Generalist"
		# 	),
		# 	vecpar = TRUE,
		# 	eval.only = TRUE,
		# 	control = mle2.control #,
		# 	# ...
		# )

		# if(modeltype == "Holling II Generalist Generalist"){
		# 	print(hollingII.via.mle2@coef)
		# 	return(hollingII.via.mle2)
		# }else{
		# 	start <- list(
		# 		attack_i = coef(hollingII.via.mle2)["attack_i"],
		# 		attack_j = coef(hollingII.via.mle2)["attack_j"],
		# 		handling_i = coef(hollingII.via.mle2)["handling_i"],
		# 		handling_j = coef(hollingII.via.mle2)["handling_j"]
		# 	)

			if(modeltype == "Holling II Hybrid Specialist" | modeltype == "Holling II Hybrid Generalist" | modeltype == "Holling II Hybrid Hybrid"){
				start$phi_ij <- 0.5
				if(phi.transform(1)==exp(1)) start$phi_ij <- log(start$phi_ij)
			}

			if(modeltype == "Holling II Specialist Hybrid" | modeltype == "Holling II Generalist Hybrid" | modeltype == "Holling II Hybrid Hybrid"){
				start$phi_ji <- 0.5
				if(phi.transform(1)==exp(1)) start$phi_ji <- log(start$phi_ji)
			}

			# fit more complex Holling Type model with above Type II as starting parameter values
			fit.via.sbplx <- nloptr::sbplx(
				x0 = unlist(start),
				fn = holling.like.1pred.2prey.NLL,
				Ni = d$Nprey1,
				Nj = d$Nprey2,
				Ni_consumed = d$Nconsumed1,
				Nj_consumed = d$Nconsumed2,
				Npredators = d$Npredator,
				time=d$Time,
				replacement = s$replacement,
				modeltype = modeltype,
				phi.transform = phi.transform,
				control = nloptr.control #,
				# ...
			)

			mle2.start <- as.list(fit.via.sbplx$par)
			names(mle2.start) <- names(start)

			# refit this model with mle2
			fit.via.mle2 <- bbmle::mle2(
				holling.like.1pred.2prey.NLL,
				start=mle2.start,
				data=list(
					Ni = d$Nprey1,
					Nj = d$Nprey2,
					Ni_consumed = d$Nconsumed1,
					Nj_consumed = d$Nconsumed2,
					Npredators = d$Npredator,
					time=d$Time,
					replacement = s$replacement,
					modeltype = modeltype,
					phi.transform = phi.transform
				),
				vecpar = TRUE,
				# eval.only = TRUE,
				control = mle2.control #,
				# ...
			)

			fit.via.mle2@call$data$sbplx.start <- unlist(start)
			names(fit.via.mle2@call$data$sbplx.start) <- names(start)

			# print(fit.via.mle2@coef)
			return(fit.via.mle2)

			# DEBUG NONE OF THE BELOW IS CURRENTLY BEING RUN!

			# convert mle2 estimation to list of starting values
			mle2.start <- as.list(fit.via.mle2@coef)
			names(mle2.start) <- names(start)

			# apparently this helps optimize over complex likelihood surfaces and get SEs when they weren't there otherwise...
			mle2.control$parscale <- abs(fit.via.mle2@coef)

			# refit with mle2 using parscale to help get an appropriate covariance matrix for the parameters
			refit.via.mle2 <- bbmle::mle2(
				holling.like.1pred.2prey.NLL,
				start=mle2.start,
				data=list(
					Ni = d$Nprey1,
					Nj = d$Nprey2,
					Ni_consumed = d$Nconsumed1,
					Nj_consumed = d$Nconsumed2,
					Npredators = d$Npredator,
					time=d$Time,
					replacement = s$replacement,
					modeltype = modeltype,
					phi.transform = phi.transform
				),
				vecpar = TRUE,
				control = mle2.control				
			)

			return(refit.via.mle2)
		# }
	}
}
