
# generate boostrapped data from a set of means and SEs of experimental observations
bootstrap.data <- function(rawdata, expttype=c("integrated","replacement")){
	# make sure that the expttype variable is valid
	expttype <- match.arg(expttype)

	# determine the "names" of all prey present
	prey <- grep("[.]mean$", colnames(rawdata), value=TRUE)
	prey <- gsub("[.]mean$", "", prey)
	prey <- gsub("^Nconsumed", "", prey)

	# generate a container for the data
	if("Time" %in% colnames(rawdata)){
		d <- matrix(NA, 0, 2 + 2 * length(prey))
	}else{
		d <- matrix(NA, 0, 1 + 2 * length(prey))
	}

	# go row by row through the data
	for(r in 1:nrow(rawdata)){
		# if for some reason we do not know the SE for any of the different observations treat whole set as a SINGLE replicate
		if(any(is.na(rawdata[r,paste0("Nconsumed",prey,".se")]))){
			rawdata[r,paste0("Nconsumed",prey,".se")] <- 0
			rawdata[r,"n"] <- 1
		}

		# go replicate by replicate
		for(e in 1:rawdata[r,"n"] ){
			Na <- c()
			Ne <- c()
			# sample prey by prey
			for(i in prey){
				# known number of available prey
				Na <- c(
					Na, rawdata[r,paste0("Nprey",i)]
				)

				# sample number of eaten prey
				Ne <- c(
					Ne,
					fr.sample(
						mean=rawdata[r,paste0("Nconsumed",i,".mean")],
						se=rawdata[r,paste0("Nconsumed",i,".se")],
						n=rawdata[r,"n"],
						expttype=expttype,
						Nprey=rawdata[r,paste0("Nprey",i)]
					)
				)
			}

			# add the requisite infoto the new row
			new.row <- c(
				rawdata[r,"Npredator"],
				Na,
				Ne
			)
			
			# if there was time/duration data in the original data frame, put it back in
			if("Time" %in% colnames(rawdata)){
				new.row <- c(
					new.row,
					rawdata[r,"Time"]
				)
			}

			# add the new "observation" to the data set
			d <- rbind(d, new.row)
		}
	}

	# all data frames should emerge from this function with identical column names
	d <- data.frame(d, row.names=seq(1,nrow(d)))
	if("Time" %in% colnames(rawdata)){
		colnames(d) <- c(
			"Npredator",
			paste0("Nprey", prey),
			paste0("Nconsumed", prey),
			"Time"
		)
	}else{
		colnames(d) <- c(
			"Npredator",
			paste0("Nprey", prey),
			paste0("Nconsumed", prey)
		)
	}

	return(d)
}

# sample a number consumed given mean, se, n, and Nprey (which is necessary for non-replacement studies)
fr.sample <- function(mean, se, n, expttype=c("integrated", "replacement"), Nprey=NULL){
	# make sure that the expttype variable is valid
	expttype <- match.arg(expttype)

	# we require Nprey for non-replacement since they are proportion of a total which must be specified
	if(is.null(Nprey) && expttype == "integrated"){
		stop("fr.sample requires Nprey argument for non-replacement experiment types")
	}

	# normal distribution nomenclature
	mu <- mean
	sigma <- se * sqrt(n)

	if(expttype=="integrated"){
		# calculate the proportion of prey consumed; do not divide by zero prey since no prey available means no prey consumed!
		prob <- ifelse(
			Nprey==0,
			0,
			rnorm(n=1, mean=mu, sd=sigma) / Nprey
		)

		# make sure we end up with a valid proportion
		prob <- max(0, min(1, prob))
		
		# bootstrapped number eaten is the result of a binomial process
		Ne <- rbinom(
				n=1,
				size=Nprey,
				p=prob
		)
	}

	if(expttype=="replacement"){
		# calculate the expected number consumed given the rate of consumption
		lambda <- rnorm(n=1, mean=mu, sd=sigma)

		# make sure we end up with a valid rate (>=0)
		lambda <- max(0, lambda)

		# boostrapped number eaten is the result of a poisson process
		Ne <- rpois(n=1, lambda=lambda)
	}

	return(Ne)
}
