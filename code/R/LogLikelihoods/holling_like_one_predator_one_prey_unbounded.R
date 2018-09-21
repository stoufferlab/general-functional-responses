
#############################################
# holling-like functional responses
#############################################

# libraries required below
library(lamW)

# predicted number of species consumed given parameters of a holling-like functional response
holling.like.1pred.1prey = function(N0, a, h, c, phi_numer, phi_denom, P, T, expttype=c("integrated","replacement"), Pminus1=c(TRUE,FALSE)){
	expttype <- match.arg(expttype)

	# if only P-1 individuals interference with predators that are doing the feeding
	if(Pminus1){
		P_interfering <- P - 1
	}else{
		P_interfering <- P
	}

	if(expttype=="integrated"){
		if(h==0){
			N <- N0 * (1 - exp(-a * P * T))
		}else{
			heff <- h * (1 + (1 - phi_numer * phi_denom) * c * P_interfering)
			Q <- (1 + c * P_interfering)
			X <- (1 + (1 - phi_numer) * c * P_interfering)
			N <- N0 - (Q / (a * heff)) * lamW::lambertW0(((a * heff * N0)/ Q) * exp(- (a / Q) * (X * P * T - heff * N0)))
			
			# sometimes the argument in the exponential passed to lambertW0 causes it to blow up
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

	# in a world with replacement everything is hunky dory
	if(expttype=="replacement"){
		numer <- (a * N0 * (1 + (1 - phi_numer) * c * P_interfering))
		denom <- (1 + a * h * N0 + c * P_interfering + (1 - phi_numer * phi_denom) * c * P_interfering * a * h * N0)
		N <- (numer / denom) * P * T
	}

	return(N)
}

# convenience function for these holling-like models that will allow the use of nloptr instead of mle2
holling.like.1pred.1prey.NLL = function(params, modeltype, initial, killed, predators, expttype, Pminus1, time=NULL){
	if(modeltype == "Holling I"){
		attack <- exp(params[1])
		handling <- 0
		interference <- 0
		phi_numer <- 1
		phi_denom <- 1
	}

	if(modeltype == "Holling II"){
		attack <- exp(params[1])
		handling <- exp(params[2])
		interference <- 0
		phi_numer <- 1
		phi_denom <- 1
	}

	if(modeltype == "Beddington-DeAngelis"){
		attack <- exp(params[1])
		handling <- exp(params[2])
		interference <- exp(params[3])
		phi_numer <- 1
		phi_denom <- 1
	}

	if(modeltype == "Crowley-Martin"){
		attack <- exp(params[1])
		handling <- exp(params[2])
		interference <- exp(params[3])
		phi_numer <- 1
		phi_denom <- 0
	}

	if(modeltype == "Stouffer-Novak I"){
		attack <- exp(params[1])
		handling <- exp(params[2])
		interference <- exp(params[3])
		phi_numer <- 1
		phi_denom <- params[4]
	}

	if(modeltype == "Stouffer-Novak II"){
		attack <- exp(params[1])
		handling <- exp(params[2])
		interference <- exp(params[3])
		phi_numer <- params[4]
		phi_denom <- 1
	}
	
	if(modeltype == "Stouffer-Novak III"){
		attack <- exp(params[1])
		handling <- exp(params[2])
		interference <- exp(params[3])
		phi_numer <- params[4]
		phi_denom <- params[5]
	}

	# if no times are specified then normalize to time=1
	if(is.null(time)){
		time <- 1
	}

	# expected number consumed given data and parameters
	Nconsumed <- holling.like.1pred.1prey(N0=initial, a=attack, h=handling, c=interference, phi_numer=phi_numer, phi_denom=phi_denom, P=predators, T=time, expttype=expttype, Pminus1=Pminus1)

	# DEBUG if the parameters are not biologically plausible, neither should be the likelihood
	if(any(Nconsumed <= 0) | any(is.nan(Nconsumed))){
		nll <- Inf
	}else{
		# negative log likelihood based on proportion consumed (no replacement)
		if(expttype=="integrated"){
			nll <- -sum(dbinom(killed, prob=Nconsumed/initial, size=initial, log=TRUE))
		}

		# negative log likelihood based on total number consumed (replacement)
		if(expttype=="replacement"){
			nll <- -sum(dpois(killed, Nconsumed, log=TRUE))
		}
	}

	return(nll)
}

# DEBUG there must be a more elegant way to do this "within" the function itself
parnames(holling.like.1pred.1prey.NLL) <- c(
	'attack',
	'handling',
	'interference',
	'phi_numer',
	'phi_denom'
)
