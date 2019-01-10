source('../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../lib/holling_method_one_predator_one_prey.R')
source('../lib/ratio_method_one_predator_one_prey.R')
load('../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

# ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]


fit.order <- order.of.fits(ffr.fits, model="Arditi.Akcakaya", parm="exponent" , order=TRUE)
fit.order <- 4

plot.coefs(
   ffr.fits[fit.order],
   model="Arditi.Akcakaya",
   parameter="exponent",
   ilink=identity,
   plot.SEs=TRUE,
   display.outlier.ests=TRUE,
   xlab="Arditi-Akcakaya interference rate (m)",
   ylab="Dataset",
   xlim=c(0,5)
)

# dev.copy2pdf(file="../../../results/R/OnePredOnePrey_m.pdf")
# dev.off()
