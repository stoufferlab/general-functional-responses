
#############################################
# holling-like functional responses
#############################################

# libraries used below
library(deSolve)
library(doParallel)
registerDoParallel(cores=1)

# predict the number of prey consumed in a replacement or integrated experiment
# for a single predator with multiple prey species and potential "handling-time exchangability" between those prey
holling.like.1pred.2prey = function(Ni, ai, hi, Nj, aj, hj, phi_ij, phi_ji, P, T, expttype=c("integrated","replacement")){
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
			return(N0 - Nf)
		}
	}

	return(Ne)
}

holling.like.1pred.2prey.NLL = function(params, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time=rep(1,length(Ni))){
	attack_i <- params[1]
	handling_i <- params[2]
	attack_j <- params[3]
	handling_j <- params[4]
	phi_ij <- params[5]
	phi_ji <- params[6]

	# we use parameter transformations to help improve the fitting and to avoid needing bounded optimization
	# standard Holling Type II parameters
	attack_i <- exp(attack_i)
	attack_j <- exp(attack_j)
	handling_i <- exp(handling_i)
	handling_j <- exp(handling_j)

	# new parameters
	# phi_ij <- plogis(phi_ij)
	# phi_ji <- plogis(phi_ji)

	# expected number consumed
	Nconsumed <- holling.like.1pred.2prey(
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

	# disallow biologically implausible predictions
	if(any(Nconsumed < 0 | is.nan(Nconsumed))){
		return(Inf)
	}

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

holling.like.1pred.2prey.NLL.same_phi = function(attack_i, handling_i, attack_j, handling_j, phi, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time=rep(1,length(Ni))){
	nll <- holling.like.1pred.2prey.NLL(attack_i, handling_i, attack_j, handling_j, phi, phi, Ni, Nj, Ni_consumed, Nj_consumed, Npredators, expttype, time)
	return(nll)
}
