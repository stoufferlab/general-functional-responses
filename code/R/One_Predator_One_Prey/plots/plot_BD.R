source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now profile the fits
ffr.fits <- profile_coefs(ffr.fits, 
                          model='Beddington.DeAngelis',
                          point.est='median',
                          printWarnings = TRUE,
                          which.pars = "interference")

save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.BD.Rdata')
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.BD.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General data and plot preparations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

###################################################
# ~~~~~~~~~~~~~~~~ BD interference ~~~~~~~~~~~~~~~~
###################################################

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_BD_int.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits,
     model="Beddington.DeAngelis",
     parameter="interference",
     ilink=exp,
     point.est='median',
     plot.SEs=TRUE,
     display.outlier.ests=TRUE,
     xlab="Beddington-DeAngelis interference rate (c)",
     labels=labels,
     vertLines=NA,
     xlim=c(0,8)
  )
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
int <- unlist(lapply(ffr.fits, function(x) x$estimates[['Beddington.DeAngelis']]["50%",'interference',"estimate"]))
parm <- exp(int)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_BD_int_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  ylim <- c(0,5.5)
  plot(parm~sample.sizes,
       ylim=ylim,
       type='n',
       log='x', 
       xlab='Sample size (n)', 
       ylab='Beddington-DeAngelis interference rate (c)')
  arrows(sample.sizes[parm>ylim[2]], 1*ylim[2], sample.sizes[parm>ylim[2]], 1.03*ylim[2], length=0.02)
  text(sample.sizes[parm>ylim[2]],0.98*ylim[2],round(parm[parm>ylim[2]],0), cex=0.5)
  points(sample.sizes,parm,
         pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_BD_int_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm[parm<8],breaks=50,
       xlab='Beddington-DeAngelis interference rate (c)',
       main='',
       col='grey',
       ylim=c(0,16),
       xlim=c(0,8.5))
  abline(v=1,lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()

# ~~~~~~~~~~~~~~~~~
a <- unlist(lapply(ffr.fits, function(x) x$estimates[['Beddington.DeAngelis']]["50%",'attack',"estimate"]))
h <- unlist(lapply(ffr.fits, function(x) x$estimates[['Beddington.DeAngelis']]["50%",'handling',"estimate"]))
parm <- exp(int)/(exp(a)*exp(h))

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_BD_int-ah_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  ylim <- c(1E-1,1E5)
  plot(parm~sample.sizes,
       ylim=ylim,
       type='n',
       log='xy', 
       xlab='Sample size (n)', 
       ylab='Ratio of interference rate to attack*handling')
  arrows(sample.sizes[parm>ylim[2]], 1*ylim[2], sample.sizes[parm>ylim[2]], 1.05*ylim[2], length=0.02)
  arrows(sample.sizes[parm<ylim[1]], 1*ylim[1], sample.sizes[parm<ylim[1]], 0.95*ylim[1], length=0.02)
  abline(h=1, col='grey', lty=2)
  points(sample.sizes,parm,
         pch=21, bg='grey')
dev.off()



