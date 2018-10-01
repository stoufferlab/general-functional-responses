
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

fit.holling.like <- function(d, s, modeltype, nloptr.control=list(), mle2.control=list(), ...){
	# estimate starting value from the data using linear regression
	x0 <- log(coef(lm(d$Nconsumed~0+I(d$Npredator * d$Nprey))))

	# fit Holling Type I via MLE with above starting parameter value
	hollingI.via.sbplx <- nloptr::sbplx(
		x0 = x0,
		fn = holling.like.1pred.1prey.NLL,
		modeltype="Holling I",
		initial=d$Nprey,
		killed=d$Nconsumed,
		predators=d$Npredator,
		time=d$Time,
		expttype=s$expttype,
		Pminus1=s$Pminus1,
		control = nloptr.control,
		...
	)

	# refit with mle2 since this also estimates the covariance matrix for the parameters
	hollingI.via.mle2 <- bbmle::mle2(
		minuslogl = holling.like.1pred.1prey.NLL,
		start = list(
			attack = hollingI.via.sbplx$par[1]
		),
		data = list(
			modeltype="Holling I",
			initial=d$Nprey,
			killed=d$Nconsumed,
			predators=d$Npredator,
			time=d$Time,
			expttype=s$expttype,
			Pminus1=s$Pminus1
		),
		vecpar = TRUE,
		control = mle2.control,
		...
	)

	if(modeltype == "Holling I"){
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
			modeltype = "Holling II",
			initial = d$Nprey,
			killed = d$Nconsumed,
			predators = d$Npredator,
			time = d$Time,
			expttype = s$expttype,
			Pminus1 = s$Pminus1,
			control = nloptr.control,
			...
		)

		mle2.start <- as.list(hollingII.via.sbplx$par)
		names(mle2.start) <- names(start)

		# refit with mle2 since this also estimates the covariance matrix for the parameters
		hollingII.via.mle2 <- bbmle::mle2(
			minuslogl = holling.like.1pred.1prey.NLL,
			start = mle2.start,
			data = list(
				modeltype = "Holling II",
				initial = d$Nprey,
				killed = d$Nconsumed,
				predators = d$Npredator,
				time = d$Time,
				expttype = s$expttype,
				Pminus1 = s$Pminus1
			),
			vecpar = TRUE,
			control = mle2.control
		)

		if(modeltype == "Holling II"){
			return(hollingII.via.mle2)
		}else{
			if(modeltype == "Beddington-DeAngelis"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1)
				)
			}

			if(modeltype == "Crowley-Martin"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1)
				)
			}

			if(modeltype == "Stouffer-Novak I"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1),
					phi_denom = 1
				)
			}

			if(modeltype == "Stouffer-Novak II"){
				start <- list(
					attack = coef(hollingII.via.mle2)["attack"],
					handling = log(1), #coef(fit.via.mle2)["handling"],
					interference = log(1),
					phi_numer = 1
				)
			}

			if(modeltype == "Stouffer-Novak III"){
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
				expttype = s$expttype,
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
					expttype = s$expttype,
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
					expttype = s$expttype,
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
