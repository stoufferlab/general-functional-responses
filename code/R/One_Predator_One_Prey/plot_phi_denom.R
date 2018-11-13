source('../lib/plot_coefs.R')

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

# generate a quick and dirty plot of the phi_denom parameters of the SNI model
dev.new()

plot.coefs(ffr.fits[how.to.order], #[36:40],
	"Stouffer.Novak.I",
	"phi_denom",
	plot.SEs=TRUE,
	ilink=identity,
	xlab="Effect of feeding on interfering",
	ylab="Dataset",
	xlim=c(-10,10)
)

dev.copy2pdf(file="../../../results/R/OnePredOnePrey_PhiDenom.pdf")
dev.off()
