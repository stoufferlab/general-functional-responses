
# DEBUG needs a concrete SE estimate and actual sample sizes
bootstrap.data <- function(rawdata, expttype=c("integrated","replacement")){
	expttype <- match.arg(expttype)

	# generate a container for the data
	d <- matrix(NA, sum(rawdata$n), 3)

	# bootstrap
	counter <- 1	
	for(r in 1:nrow(rawdata)){
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

			d[counter,] <- c(rawdata[r,"Npredator"],rawdata[r,"Nprey"],Ne)
			counter <- counter + 1
		}
	}
	
	d <- data.frame(d)
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed")

	return(d)
}