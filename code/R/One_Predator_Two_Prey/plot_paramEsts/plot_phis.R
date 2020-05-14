
# source('plot_coefs.R')


plot.coefs <- function(ffr.fits, modeltype, parameters, plot.SEs=FALSE, ffr.cfs=NULL, ilink=identity, xlim=NULL, ... ){
	# scaffold of a plot that doesn't actually show anything
	plot(
		y = 1:(2*length(ffr.fits)),
		x = unlist(lapply(ffr.fits, function(x){ilink(coef(x$fits[[modeltype]])[parameters])})),
		type='n',
		yaxt='n',
		xlim=xlim,
		...
	)

	# tick marks to indicate different data sets
	axis(side=2, at=1:(2*length(ffr.fits)), labels=FALSE)

	# mark where the existing models fall
	abline(v=0, lty=2)
	abline(v=1, lty=2)

	# length of arrows to indicate values off the plot
	delta.arrow <- 0.15

	# plot these bad boys
	i <- 1
	for(nn in names(ffr.fits)){
		x <- ffr.fits[[nn]]

		# color points depending on predator/parasitoid
		# col <- ifelse(x$study.info$predator, "black", "red")	
		col <- "black"

		# the median estimate is easy to determine regardless of the type of data
		mm <- x$estimates[[modeltype]]["50%",parameters,"estimate"]

		# cheeky upper and lower bounds in the absence of SE information
		lb <- mm
		ub <- mm

		# make all lines the equivalent for now
		lty <- "solid"

		# different ways to estimate the SEs from the model(s)
		if(plot.SEs){
			# if we did not need to bootstrap the data				
			if(x$estimates[[modeltype]]["n",1,1] == 1){
				# we try to precalculate the profile confidence intervals since it takes a very long time
				if(!is.null(ffr.cfs) || !is.null(ffr.cfs[[nn]][[modeltype]])){
					cf <- ffr.cfs[[nn]][[modeltype]]
				}else{
					# cf <- try(confint(x$fits[[modeltype]], try_harder=TRUE, level=0.68, tol.newmin=Inf, quietly=TRUE))
					cf <- 1
					class(cf) <- "try-error"
				}
				
				# seems like the profiling code was successful
				if(!inherits(cf, "try-error")){
					# best case equals solid line
					lty <- "solid"

					lb <- cf[parameters,1]
					ub <- cf[parameters,2]

					# sometimes we profile things but still get NA intervals
					lb <- ifelse(is.na(lb), xlim[1], lb)
					ub <- ifelse(is.na(ub), xlim[2], ub)
				}else{
					# use the quadratic assumption?

					# quadratic approximation equals dashed line
					lty <- "dashed"

					# get the SEs directly from the model output
					se <- coef(summary(x$fits[[modeltype]]))[parameters,"Std. Error"]
					lb <- ub <- se
					for(j in c(1,2)){
						lb[j] <- ifelse(is.na(se[j]), xlim[1], mm[j] - se[j])
						ub[j] <- ifelse(is.na(se[j]), xlim[2], mm[j] + se[j])
					}
				}
			}else{
				# we bootstrapped the data and will use the quantiles of those fits to match +/- 1 SE
				
				# bootstrapped is always dotted line
				# dotted
				lty <- "dotted"

				# use the central interval equivalent to one SD as the bounds
				lb <- x$estimates[[modeltype]]["16%",parameters,"estimate"]
				ub <- x$estimates[[modeltype]]["84%",parameters,"estimate"]
			}

			# DEBUG
			# # don't plot off the figure
			# lb <- ifelse(lb < xlim[1], xlim[1], lb)
			# ub <- ifelse(ub > xlim[2], xlim[2], ub)
			for(j in c(1,2)){
				lb[j] <- ifelse(lb[j] < xlim[1], xlim[1], lb[j])
				ub[j] <- ifelse(ub[j] > xlim[2], xlim[2], ub[j])
			}
		}

		print(lb)
		print(ub)
		for(j in c(1,2)){
			if(mm[j] > xlim[1] & mm[j] < xlim[2]){
				# draw the error bars
				segments(ilink(lb[j]), i, ilink(ub[j]), i, col=col, lty=lty)

				# arrow up the limiting cases
				if(lb[j]==xlim[1]){
					arrows(xlim[1], i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
				}
				if(ub[j]==xlim[2]){
					arrows(xlim[2], i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
				}
				
				# plot the actual mean estimate
				points(y=c(i),x=c(ilink(mm[j])),col=col,bg=col,pch=21)
			}else{
				if(mm[j] > xlim[2]){
					arrows(xlim[2]-delta.arrow, i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
				}else{
					arrows(xlim[1]+delta.arrow, i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
				}
			}

			i <- i + 1
		}
	}
}

xlim <- c(-5, 5)

modeltype <- "Holling II Hybrid Hybrid"
parameters <- c("phi_ij","phi_ji")

# dev.new()

plot.coefs(
	ffr.fits,
	modeltype = modeltype,
	parameters = parameters,
	plot.SEs = TRUE,
	xlim = xlim,
	ilink=identity,
	xlab="Effect of feeding on feeding",
	ylab="Dataset"
)

# dev.copy2pdf(file="../../../results/R/OnePredTwoPrey_Phis.pdf")
# dev.off()
