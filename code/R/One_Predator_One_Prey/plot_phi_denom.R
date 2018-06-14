
# source('plot_coefs.R')


plot.coefs <- function(ffr.fits, modeltype, parameter, plot.SEs=FALSE, ilink=identity, xlim=NULL, ... ){
	plot(
		y = 1:length(ffr.fits),
		x = unlist(lapply(ffr.fits, function(x){ilink(coef(x[[modeltype]])[parameter])})),
		type='n',
		yaxt='n',
		xlim=xlim,
		...
	)
	axis(side=2, at=1:length(ffr.fits), labels=FALSE)

	i <- 1
	for(nn in names(ffr.fits)){
		x <- ffr.fits[[nn]]
		mm <- max(xlim[1],min(xlim[2],coef(x[[modeltype]])[parameter]))
		mm <- coef(x[[modeltype]])[parameter]
		if(plot.SEs){
			# # xx <- try({
			# 	p0 <- profile(x[[modeltype]], try_harder=TRUE)
			# 	cnfint <- confint(p0, level=25, try_harder=TRUE)
			# 	# print(class(cnfint))
			# 	lb <- cnfint[parameter,1]
			# 	ub <- cnfint[parameter,2]
			# # })
			# xx <- 0

			# if(inherits(xx, "try-error")){
				se <- coef(summary(x[[modeltype]]))[parameter,"Std. Error"]
				if(!is.nan(se)){
					lb <- mm - se
					ub <- mm + se
				}else{
					lb <- mm
					ub <- mm
				}
			# }
			# print(se)
			# if(!is.nan(se)){
			segments(ilink(lb), i, ilink(ub), i)
			# }
		}
		points(y=c(i),x=c(ilink(mm)))
		i <- i + 1
	}

	abline(v=0)
	abline(v=1)
}


# sort by mean of phi_denom
if(TRUE){
	phi_denoms <- sapply(ffr.fits, function(x) coef(x[["Stouffer-Novak I"]])["phi_denom"])
	ffr.fits <- ffr.fits[order(phi_denoms)]
}

# ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]

# prepare some plots of the fits
# DEBUG work on the confidence intervals aspect
# par(mfrow=c(2,2))

plot.coefs(ffr.fits, #[36:40],
	"Stouffer-Novak I",
	"phi_denom",
	plot.SEs=TRUE,
	ilink=identity,
	xlab="Effect of feeding on interfering [Crowley-Martin = 0; Beddington-DeAngelis = 1]",
	ylab="Dataset",
	xlim=c(-5,5)
)
