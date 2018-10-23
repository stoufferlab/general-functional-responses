
# source('plot_coefs.R')


plot.coefs <- function(ffr.fits, modeltype, parameter, plot.SEs=FALSE, ilink=identity, xlim=NULL, ... ){
	# scaffold of a plot that doesn't actually show anything
	plot(
		y = 1:length(ffr.fits),
		x = unlist(lapply(ffr.fits, function(x){ilink(coef(x$fits[[modeltype]])[parameter])})),
		type='n',
		yaxt='n',
		xlim=xlim,
		...
	)

	# tick marks to indicate different data sets
	axis(side=2, at=1:length(ffr.fits), labels=FALSE)

	# mark where the existing models fall
	abline(v=0, lty=2)
	abline(v=1, lty=2)

	# length of arrows to indicate values off the plot
	delta.arrow <- 0.1

	# plot these bad boys
	i <- 1
	for(nn in names(ffr.fits)){
		x <- ffr.fits[[nn]]

		# color points depending on predator/parasitoid
		col <- ifelse(x$study.info$predator, "black", "red")	

		# the median estimate is easy to determine regardless of the type of data
		mm <- x$estimates[[modeltype]]["50%",parameter,"estimate"]

		# cheeky upper and lower bounds in the absence of SE information
		lb <- mm
		ub <- mm

		# make all lines the equivalent for now
		lty <- "solid"

		# different ways to estimate the SEs from the model(s)
		if(plot.SEs){
			# if we did not need to bootstrap the data				
			if(x$estimates[[modeltype]]["n",1,1] == 1){
				# estimate the profile confidence interval in first instance
				cf <- try(confint(x$fits[[modeltype]], try_harder=TRUE, level=0.68, tol.newmin=Inf, quietly=TRUE))
				
				# seems like the profiling code was successful
				if(!inherits(cf, "try-error")){
					# best case equals solid line
					lty <- "solid"

					lb <- cf[parameter,1]
					ub <- cf[parameter,2]

					# sometimes we profile things but still get NA intervals
					lb <- ifelse(is.na(lb), xlim[1], lb)
					ub <- ifelse(is.na(ub), xlim[2], ub)
				}else{
					# use the quadratic assumption?

					# quadratic approximation equals dashed line
					lty <- "dashed"

					# get the SEs directly from the model output
					se <- coef(summary(x$fits[[modeltype]]))[parameter,"Std. Error"]
					lb <- ifelse(is.na(se), xlim[1], mm - se)
					ub <- ifelse(is.na(se), xlim[2], mm + se)
				}
			}else{
				# we bootstrapped the data and will use the quantiles of those fits to match +/- 1 SE
				
				# bootstrapped is always dotted line
				# dotted
				lty <- "dotted"

				# use the central interval equivalent to one SD as the bounds
				lb <- x$estimates[[modeltype]]["16%",parameter,"estimate"]
				ub <- x$estimates[[modeltype]]["84%",parameter,"estimate"]
			}

			# don't plot off the figure
			lb <- ifelse(lb < xlim[1], xlim[1], lb)
			ub <- ifelse(ub > xlim[2], xlim[2], ub)
		}

		if(mm > xlim[1] & mm < xlim[2]){
			# draw the error bars
			segments(ilink(lb), i, ilink(ub), i, col=col, lty=lty)

			# arrow up the limiting cases
			if(lb==xlim[1]){
				arrows(xlim[1], i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
			}
			if(ub==xlim[2]){
				arrows(xlim[2], i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
			}
			
			# plot the actual mean estimate
			points(y=c(i),x=c(ilink(mm)),col=col,bg=col,pch=21)
		}else{
			if(mm > xlim[2]){
				arrows(xlim[2]-delta.arrow, i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
			}else{
				arrows(xlim[1]+delta.arrow, i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
			}
		}

		i <- i + 1
	}
}


# sort by mean of phi_denom
if(FALSE){
	phi_denoms <- unlist(lapply(ffr.fits, function(x) coef(x$fits[["Stouffer.Novak.I"]])["phi_denom"]))
	how.to.order <- order(phi_denoms)
	# ffr.fits <- ffr.fits[names(phi_denoms)[order(phi_denoms)]]
}else{
	how.to.order <- 1:length(ffr.fits)
}

# ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]

# prepare some plots of the fits
# DEBUG work on the confidence intervals aspect
# par(mfrow=c(2,2))

plot.coefs(ffr.fits[how.to.order], #[36:40],
	"Stouffer.Novak.I",
	"phi_denom",
	plot.SEs=TRUE,
	ilink=identity,
	xlab="Effect of feeding on interfering",
	ylab="Dataset",
	xlim=c(-2,2)
)
