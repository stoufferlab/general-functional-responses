source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ffr.fits <- profile_coefs(ffr.fits, 
                          model='Stouffer.Novak.I',
                          point.est='median',
                          printWarnings = TRUE)

save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General data and plot preparations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Stouffer.Novak.I", order.parm="phi_denom")
ffr.fits <- ffr.fits[fit.order]

# fraction of replicates in which significant depletion occurred (in non-replacement datasets)
col.vec<-rep('black',length(ffr.fits))
depleted <- unlist(lapply(ffr.fits, depletion.check, cutoff=0.7))
col.vec[depleted>0] <- 'red'

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

###################################################
# ~~~~~~~~~~~~~~~~~~ SN1 PhiDenom ~~~~~~~~~~~~~~~~~
###################################################
pdf(file="../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_SN1_PhiDenom.pdf",height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
	model="Stouffer.Novak.I",
	parameter="phi_denom",
  ilink=identity,
  point.est='median',
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
	xlab="Effect of feeding on interfering",
  labels=labels,
  vertLines = c(0,1),
	xlim=c(-3,3)
)
dev.off()

