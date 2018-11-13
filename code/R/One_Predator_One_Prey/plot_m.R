source('../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
load('../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

# ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]

# prepare some plots of the fits
# DEBUG work on the confidence intervals aspect
# par(mfrow=c(2,2))

# generate a quick and dirty plot of the phi_denom parameters of the SNI model
dev.new()

fit.order <- order.of.fits(ffr.fits, model="Arditi-Akcakaya", parm="m" , order=FALSE)

plot.coefs(ffr.fits[fit.order], #[36:40],
   "Arditi-Akcakaya",
   "m",
   plot.SEs=TRUE,
   ilink=identity,
   xlab="Interference rate",
   ylab="Dataset",
   xlim=c(0,5)
)

dev.copy2pdf(file="../../../results/R/OnePredOnePrey_m.pdf")
dev.off()
