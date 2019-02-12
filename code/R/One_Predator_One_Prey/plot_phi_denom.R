source('../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../lib/holling_method_one_predator_one_prey.R')
source('../lib/ratio_method_one_predator_one_prey.R')

load('../../../results/R/ffr.fits_OnePredOnePrey.Rdata')


fit.order <- order.of.fits(ffr.fits, model="Stouffer.Novak.I", order.parm="phi_denom" , order=TRUE)
ffr.fits<-ffr.fits[fit.order]

pdf(file="../../../results/R/OnePredOnePrey_PhiDenom.pdf",height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
	"Stouffer.Novak.I",
	"phi_denom",
	plot.SEs=TRUE,
	ilink=identity,
  point.est='median',
	xlab="Effect of feeding on interfering",
	ylab="",
	xlim=c(-10,10)
)
dev.off()
