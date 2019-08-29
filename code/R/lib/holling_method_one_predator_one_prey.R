#############################################
# holling-like functional responses
#############################################
# For non-replacement datasets, we have two options : 
# (1) solve using lambertsW (or solve transcendental eqn directly if needed)
# (2) integrate
#############################################
# libraries required
library(bbmle)
library(nloptr)
library(lamW)
library(odeintr)

sp <- list.files("../../..", "set_params.R", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
source(sp)
#############################################

# For integration method, define the ode in C++ format
# holling.like.1pred.1prey.sys = '
#   // generalized functional response for one predator one prey
#     dxdt[0] = -P * (a * x[0] * (1 + (1 - phi_numer) * c * P_interfering)) / (1 + a * h * x[0] + c * P_interfering + (1 - phi_numer * phi_denom) * c * P_interfering * a * h * x[0]);
#   
#   // consumption rate cannot be positive
#     if(dxdt[0] > 0) dxdt[0] = 0;
# '
# 
# # compile the above C++ code into something we can run in R
# odeintr::compile_sys(
#   "hl_1pred_1prey",
#   holling.like.1pred.1prey.sys,
#   pars = c("a", "h", "c", "phi_numer", "phi_denom", "P", "P_interfering") #,
#   # method = "bsd"
# )


# predicted number of species consumed given parameters of a holling-like functional response
holling.like.1pred.1prey = function(N0, a, h, c, phi_numer, phi_denom, P, T, 
                                    replacement, 
                                    Pminus1=c(TRUE,FALSE),
                                    integrate=FALSE,
                                    overrideTranscendental=FALSE){

	# if only P-1 individuals interference with predators that are doing the feeding
	if(Pminus1){
		P_interfering <- P - 1
	}else{
		P_interfering <- P
	}

	# in a world with replacement everything is hunky dory
	if(replacement){
		numer <- (a * N0 * (1 + (1 - phi_numer) * c * P_interfering))
		denom <- (1 + a * h * N0 + c * P_interfering + (1 - phi_numer * phi_denom) * c * P_interfering * a * h * N0)
		N <- (numer / denom) * P * T
		N <- pmax(0,N)
		return(N)
	}

	# without replacement
	if(!replacement){
		if(h==0){ # For Type I things are simple:
			N <- N0 * (1 - exp(-a * P * T))
		}else{ # For all other models...
  		  if(integrate){  # solve by direct integration
  		    N <- numeric(length(N0))
  		    for(i in seq.int(length(N0))){
  
  		      # set parameters within ode solver
  		      hl_1pred_1prey_set_params(a=a, h=h, c=c, 
  		                                phi_numer=phi_numer, 
  		                                phi_denom=phi_denom, 
  		                                P=P[i], 
  		                                P_interfering=P_interfering[i])
  		      
  		      # calculate the final number of prey integrating the ode
  		      Nfinal <- hl_1pred_1prey(N0[i], T[i], T[i]/1000.)
  		      
  		      # we only need the last row since this is the final "abundance"
  		      Nfinal <- as.numeric(Nfinal[nrow(Nfinal),2])

  		      # the number consumed is the difference between what we started with and what is left
  		      N[i] <- N0[i] - Nfinal
  		    }
  		  } else {	# solve using lambertsW (or trancendental equation)
      			heff <- h * (1 + (1 - phi_numer * phi_denom) * c * P_interfering)
      			Q <- (1 + c * P_interfering)
      			X <- (1 + (1 - phi_numer) * c * P_interfering)
      			N <- N0 - (Q / (a * heff)) * lamW::lambertW0(((a * heff * N0)/ Q) * exp(- (a / Q) * (X * P * T - heff * N0)))
      			
      			# sometimes the argument in the exponential passed to lambertW0 causes it to blow up
      			if(!overrideTranscendental){
      			  if(any(is.infinite(N))){
      				# the explicit result of the analytical integration without solving for N implictly
      				ffff <- function(N, N0, P, T, a, heff, Q, X){
      					dN <- Q * log((N0 - N)/N0) - a * heff * N
      					dt <- - a * X * P * T
      					dN - dt
      				}
      				# sometimes the time argument is a constant and not a vector
      				if(length(T)==1){
      					T <- rep(T, length(N0))
      				}
      				# check which predictions are non-sensical
      				for(i in 1:length(N0)){
      					if(is.infinite(N[i])){
      						# we need to solve the transcendental equation directly
      						nn <- uniroot(ffff, lower=0, upper=N0[i], N0=N0[i], P=P[i], T=T[i], a=a, heff=heff[i], Q=Q[i], X=X[i])
      						N[i] <- nn$root
      					}
      				}
      			  }
      			}
  		  }
		}
		return(N)
	}

	stop()
}

# negative log likelihood for holling-like models given parameters and requisite data
holling.like.1pred.1prey.NLL = function(params,
                                        modeltype, 
                                        initial, 
                                        killed, 
                                        predators, 
                                        replacement, 
                                        Pminus1, 
                                        time=NULL){
  set_params(params, modeltype)

	# if no times are specified then normalize to time=1
	if(is.null(time)){
		time <- 1
	}

	# expected number consumed given data and parameters
	Nconsumed <- holling.like.1pred.1prey(N0=initial,
	                                      a=attack,
	                                      h=handling,
	                                      c=interference,
	                                      phi_numer=phi_numer,
	                                      phi_denom=phi_denom,
	                                      P=predators,
	                                      T=time,
	                                      replacement=replacement,
	                                      Pminus1=Pminus1)
	
	# reduce to unique data rows to speed up. There's probably an even faster way, but...
	# d.ori <- data.frame(initial, predators, time)
	# d.uniq <- unique(d.ori)
	# Nconsumed.uniq <- holling.like.1pred.1prey(N0=d.uniq$initial,
	#                                            a=attack, 
	#                                            h=handling, 
	#                                            c=interference, 
	#                                            phi_numer=phi_numer, 
	#                                            phi_denom=phi_denom, 
	#                                            P=d.uniq$predators, 
	#                                            T=d.uniq$time, 
	#                                            replacement=replacement, 
	#                                            Pminus1=Pminus1)
	# Nconsumed <- merge(d.ori, cbind(d.uniq, Nconsumed.uniq))$Nconsumed.uniq

	# if the parameters are not biologically plausible, neither should be the likelihood
	if(any(Nconsumed <= 0) | any(is.nan(Nconsumed))){
		nll <- Inf
		return(nll)
	}else{
		# negative log likelihood based on proportion consumed (no replacement)
		if(!replacement){
		  # warnings suppressed because direct integration can return prob = 0 or 1, which results in NaNs
			nll <- suppressWarnings( -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE)) )
			if(is.nan(nll)){
			  nll <- Inf
			}
			return(nll)
		}

		# negative log likelihood based on total number consumed (replacement)
		if(replacement){
			nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
			return(nll)
		}
	}

	stop()
}

# needed by mle2 to pass named parameters in the right order
# DEBUG there must be a more elegant way to do this "within" the function itself
parnames(holling.like.1pred.1prey.NLL) <- c(
	'attack',
	'handling',
	'interference',
	'phi_numer',
	'phi_denom'
)

# given data (d), study info (s), and modeltype (e.g., "Holling I"), fit functional response data
fit.holling.like <- function(d, s, 
                             modeltype, 
                             nloptr.control=list(), 
                             mle2.control=list(), 
                             ...){

	# estimate starting value from the data using linear regression
  start <- list(
    attack = log(coef(lm(d$Nconsumed~0+I(d$Npredator * d$Nprey))))
  )

	# fit Holling Type I via MLE with above starting parameter value
	hollingI.via.sbplx <- nloptr::sbplx(
		x0 = unlist(start),
		fn = holling.like.1pred.1prey.NLL,
		modeltype="Holling.I",
		initial=d$Nprey,
		killed=d$Nconsumed,
		predators=d$Npredator,
		time=d$Time,
		replacement=s$replacement,
		Pminus1=s$Pminus1,
		control = nloptr.control,
		...
	)

	# refit with mle2 since this also estimates the covariance matrix for the parameters
	mle2.start <- as.list(hollingI.via.sbplx$par)
	names(mle2.start) <- names(start)
	
	hollingI.via.mle2 <- bbmle::mle2(
		minuslogl = holling.like.1pred.1prey.NLL,
		start = mle2.start,
		data = list(
			modeltype="Holling.I",
			initial=d$Nprey,
			killed=d$Nconsumed,
			predators=d$Npredator,
			time=d$Time,
			replacement=s$replacement,
			Pminus1=s$Pminus1
		),
		vecpar = TRUE,
		control = mle2.control,
		...
	)

	if(modeltype == "Holling.I"){
		return(hollingI.via.mle2)
	}
	# code to fit subsequent models
	else{
		# fit Holling II first
		start <- list(
			attack = coef(hollingI.via.mle2)["attack"],
			handling = log(1)
		)

		# fit the more complex model with sbplx first
		hollingII.via.sbplx <- nloptr::sbplx(
			x0 = unlist(start),
			fn = holling.like.1pred.1prey.NLL,
			modeltype = "Holling.II",
			initial = d$Nprey,
			killed = d$Nconsumed,
			predators = d$Npredator,
			time = d$Time,
			replacement = s$replacement,
			Pminus1 = s$Pminus1,
			control = nloptr.control
			# ,
			# ...
		)

		mle2.start <- as.list(hollingII.via.sbplx$par)
		names(mle2.start) <- names(start)

		# refit with mle2 since this also estimates the covariance matrix for the parameters
		hollingII.via.mle2 <- bbmle::mle2(
			minuslogl = holling.like.1pred.1prey.NLL,
			start = mle2.start,
			data = list(
				modeltype = "Holling.II",
				initial = d$Nprey,
				killed = d$Nconsumed,
				predators = d$Npredator,
				time = d$Time,
				replacement = s$replacement,
				Pminus1 = s$Pminus1
			),
			vecpar = TRUE,
			control = mle2.control
		)

		if(modeltype == "Holling.II"){
			return(hollingII.via.mle2)
		}else{
			if(modeltype == "Beddington.DeAngelis"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1)
				)
			}

			if(modeltype == "Crowley.Martin"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1)
				)
			}

			if(modeltype == "Stouffer.Novak.I"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1),
					phi_denom = 1
				)
			}

			if(modeltype == "Stouffer.Novak.II"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1),
					phi_numer = 1
				)
			}

			if(modeltype == "Stouffer.Novak.III"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1),
					phi_numer = 1,
					phi_denom = 1
				)
			}

			# fit the more complex model with sbplx first
			fit.via.sbplx <- nloptr::sbplx(
				x0 = unlist(start),
				fn = holling.like.1pred.1prey.NLL,
				modeltype = modeltype,
				initial = d$Nprey,
				killed = d$Nconsumed,
				predators = d$Npredator,
				time = d$Time,
				replacement = s$replacement,
				Pminus1 = s$Pminus1,
				control = nloptr.control,
				...
			)

			# convert nloptr estimation to list of starting values
			mle2.start <- as.list(fit.via.sbplx$par)
			names(mle2.start) <- names(start)
			
			# fit with mle2 since this provides other convenience estimates
			fit.via.mle2 <- bbmle::mle2(
				minuslogl = holling.like.1pred.1prey.NLL,
				start = mle2.start,
				data = list(
					modeltype = modeltype,
					initial = d$Nprey,
					killed = d$Nconsumed,
					predators = d$Npredator,
					time = d$Time,
					replacement = s$replacement,
					Pminus1 = s$Pminus1
				),
				vecpar = TRUE,
				control = mle2.control
			)

			# return(fit.via.mle2)

			# hess <- numDeriv::hessian(
			# 	func = holling.like.1pred.1prey.NLL,
			# 	x = fit.via.mle2@coef,
			# 	modeltype = modeltype,
			# 	initial = d$Nprey,
			# 	killed = d$Nconsumed,
			# 	predators = d$Npredator,
			# 	time = d$Time,
			# 	expttype = s$expttype,
			# 	Pminus1 = s$Pminus1
			# )
			# print(hess)

			# convert mle2 estimation to list of starting values
			mle2.start <- as.list(fit.via.mle2@coef)
			names(mle2.start) <- names(start)

			# apparently this helps optimize over complex likelihood surfaces and get SEs when they weren't there otherwise...
			mle2.control$parscale <- abs(fit.via.mle2@coef)

			# refit with mle2 using parscale to help get an appropriate covariance matrix for the parameters
			refit.via.mle2 <- bbmle::mle2(
				minuslogl = holling.like.1pred.1prey.NLL,
				start = mle2.start,
				data = list(
					modeltype = modeltype,
					initial = d$Nprey,
					killed = d$Nconsumed,
					predators = d$Npredator,
					time = d$Time,
					replacement = s$replacement,
					Pminus1 = s$Pminus1
				),
				vecpar = TRUE,
				control = mle2.control				
			)

			return(refit.via.mle2)

			# # DEBUG to solve some issues with the fitting
			# # if we couldn't get the SEs correctly, try a few more things
			# print(is.na(sqrt(diag(fit.via.mle2@vcov))))
			# if(any(is.nan(sqrt(diag(fit.via.mle2@vcov))) | is.na(sqrt(diag(fit.via.mle2@vcov))))){
			# 	# refit with mle2 since this also estimates the covariance matrix for the parameters
			# 	fit.via.mle2 <- bbmle::mle2(
			# 		minuslogl = holling.like.1pred.1prey.NLL,
			# 		start = mle2.start,
			# 		fixed = fixed,
			# 		data = list(
			# 			modeltype = modeltype,
			# 			initial = d$Nprey,
			# 			killed = d$Nconsumed,
			# 			predators = d$Npredator,
			# 			time = d$Time,
			# 			expttype = s$expttype,
			# 			Pminus1 = s$Pminus1
			# 		),
			# 		vecpar = TRUE,
			# 		control = mle2.control,
			# 		method="Nelder-Mead"
			# 	)

			# 	print(attr(fit.via.mle2, "vcov"))

			# 	hess <- numDeriv::hessian(
			# 		func = holling.like.1pred.1prey.NLL,
			# 		x = fit.via.mle2@coef,
			# 		modeltype = modeltype,
			# 		initial = d$Nprey,
			# 		killed = d$Nconsumed,
			# 		predators = d$Npredator,
			# 		time = d$Time,
			# 		expttype = s$expttype,
			# 		Pminus1 = s$Pminus1
			# 	)
			# 	print(hess)

			# 	fit.via.mle2@vcov <- MASS::ginv(hess)
			# 	# print(hess)
			# 	print(fit.via.mle2@vcov)

			# }

			
			# 	# fit the model with cobyla first
			# 	fit.via.cobyla <- nloptr::cobyla(
			# 		x0 = unlist(cobyla.start),
			# 		fn = holling.like.1pred.1prey.NLL,
			# 		# lower = lower,
			# 		# upper = upper,
			# 		modeltype = modeltype,
			# 		initial = d$Nprey,
			# 		killed = d$Nconsumed,
			# 		predators = d$Npredator,
			# 		time = d$Time,
			# 		expttype = s$expttype,
			# 		Pminus1 = s$Pminus1,
			# 		control = nloptr.control,
			# 		...
			# 	)

			
		}
	}
}
