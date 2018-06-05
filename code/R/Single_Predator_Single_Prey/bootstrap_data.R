
# DEBUG needs a concrete SE estimate and actual sample sizes
bootstrap.data <- function(rawdata, expttype=c("integrated","replacement")){
	expttype <- match.arg(expttype)

	# generate a container for the data
	if("Time" %in% colnames(rawdata)){
		d <- matrix(NA, 0, 4)
	}else{
		d <- matrix(NA, 0, 3)
	}

	

	# bootstrap
	counter <- 1	
	for(r in 1:nrow(rawdata)){
		# DEBUG: if the SE is NA then do not sample that set of replicates more than once
		if(is.na(rawdata[r,"Nconsumed.se"])){
			rawdata[r,"Nconsumed.se"] <- 0
			rawdata[r,"n"] <- 1
		}
		
		for( e in 1:rawdata[r,"n"] ){
			mu <- rawdata[r,"Nconsumed.mean"]
			sigma <- rawdata[r,"Nconsumed.se"] * sqrt(rawdata[r,"n"])

			if(expttype=="integrated"){
				prob <- rnorm(n=1, mean=mu, sd=sigma) / rawdata[r,"Nprey"]
				prob <- max(0, min(1, prob))
				Ne <- rbinom(
						n=1,
						size=rawdata[r,"Nprey"],
						p=prob
				)
			}

			if(expttype=="replacement"){
				lambda <- rnorm(n=1, mean=mu, sd=sigma)
				lambda <- max(0, lambda)
				Ne <- rpois(n=1, lambda=lambda)
			}

			if("Time" %in% colnames(rawdata)){
				d <- rbind(d, c(rawdata[r,"Npredator"],rawdata[r,"Nprey"],Ne,rawdata[r,"Time"]))
			}else{
				d <- rbind(d, c(rawdata[r,"Npredator"],rawdata[r,"Nprey"],Ne))
			}
		}
	}
	
	d <- data.frame(d)
	if("Time" %in% colnames(rawdata)){
		colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
	}else{
		colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
	}

	return(d)
}
