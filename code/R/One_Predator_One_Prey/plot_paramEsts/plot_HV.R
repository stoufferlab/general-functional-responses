source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ffr.fits <- profile_coefs(ffr.fits, 
                          model='Hassell.Varley',
                          point.est='median',
                          printWarnings = TRUE,
                          which.pars = "exponent")

save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.HV.Rdata')
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.HV.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General data and plot preparations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Hassell.Varley", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

###################################################
# ~~~~~~~~~~~~~~~~~~ HV exponent ~~~~~~~~~~~~~~~~~~
###################################################
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_HV_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
  model="Hassell.Varley",
  parameter="exponent",
  ilink=exp,
  point.est='median',
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
  color.factor='Parasitoids',
  xlab="Hassel-Varley interference rate (m)",
  labels=labels,
  xlim=c(0,5)
)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Hassell.Varley']]["50%",'exponent',"estimate"]))
parm <- exp(parm)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_HV_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  ylim <- c(0,3)
  plot(parm~sample.sizes,
       ylim=ylim,
       type='n',
       log='x', 
       xlab='Sample size (n)', 
       ylab='Hassel-Varley interference rate (m)')
  abline(h=c(0,1),lty=2,col='grey')
  points(sample.sizes,parm,
         pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_HV_m_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm,breaks=50,
       xlab='Hassel-Varley interference rate (m)',
       main='',
       col='grey',
       ylim=c(0,11),
       xlim=c(0,3))
  abline(v=1,lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()
