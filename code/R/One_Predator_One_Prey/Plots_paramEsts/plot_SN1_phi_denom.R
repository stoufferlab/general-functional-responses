source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

load('../../../../results/R/ffr.fits_OnePredOnePrey.Rdata')


fit.order <- order.of.fits(ffr.fits, model="Stouffer.Novak.I", order.parm="phi_denom", order=TRUE, point.est='median')
ffr.fits<-ffr.fits[fit.order]

# fraction of replicates in which significant depletion occurred (in non-replacement datasets)
col.vec<-rep('black',length(ffr.fits))
depleted <- unlist(lapply(ffr.fits, depletion.check, cutoff=0.7))
col.vec[depleted>0] <- 'red'

pdf(file="../../../../results/R/OnePredOnePrey_SN1_PhiDenom.pdf",height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
	model="Stouffer.Novak.I",
	parameter="phi_denom",
  ilink=identity,
  point.est='median',
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
  color.factor='None', # 'None', 'Parasitoids' or 'Replacement'
  # color.vector=col.vec, # delete or specify above plot()
	xlab="Effect of feeding on interfering",
	ylab="",
  labels=TRUE,
	xlim=c(-6,6)
)
dev.off()
