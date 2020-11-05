source('../../lib/plot_coefs.R')
source('../../lib/profile_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Remove datasets where AA2 method was NOT successfully applied
yesAA2 <- unlist(lapply(ffr.fits, function(x) 'Arditi.Akcakaya.Method.2' %in% names(x$fits)))
ffr.fits <- ffr.fits[which(yesAA2)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Profile the fits
ffr.fits <- profile_coefs(ffr.fits, 
                          model='Arditi.Akcakaya.Method.2',
                          point.est='median',
                          printWarnings = TRUE)

save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.AA2.Rdata')
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.AA2.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya.Method.2", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m2 <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya.Method.2']]["50%",'exponent',"estimate"]))
parm <- m2 
range(m2)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA2_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits,
     model="Arditi.Akcakaya.Method.2",
     parameter="exponent",
     ilink=identity,
     plot.SEs=TRUE,
     point.est='median',
     display.outlier.ests=TRUE,
     xlab="Arditi-Akcakaya interference strength (m) (Method 2)",
     labels=labels,
     xlim=c(-0.8,1.85)
  )
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Additional summary plots

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA2_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  ylim <- c(0,1.6)
  plot(parm~sample.sizes,
       ylim=ylim,
       type='n',
       log='x', 
       xlab='Sample size (n)', 
       ylab='Arditi-Akcakaya interference strength (m) (AA2)')
  abline(h=c(0,1),lty=2,col='grey')
  points(sample.sizes, parm,
         pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA2_m_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm[parm<2],
       breaks=50,
       xlab='Arditi-Akcakaya interference strength (m) (AA2)',
       main='',
       col='grey',
       ylim=c(0,3.5),
       xlim=c(-0.6,1.6))
  abline(v=c(0,1),lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()

##################################################
# Compare AA2 estimates to those of non-AA2 method
##################################################
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')

# Remove datasets where AA2 method was NOT successfully applied
yesAA2 <- unlist(lapply(ffr.fits, function(x) 'Arditi.Akcakaya.Method.2' %in% names(x$fits)))
ffr.fits <- ffr.fits[which(yesAA2)]

m <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))
m <- exp(m)
m2 <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya.Method.2']]["50%",'exponent',"estimate"]))

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA2_m_vs_AA_m.pdf',height=3,width=3)
par(pty='s',mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i', xaxs='i')
  xylim<-c(-0.5,1.5)
  plot(m,m2,
       type='n',
       xlim=xylim,
       ylim=xylim,
       xlab='Maximum likelihood estimate',
       ylab='AA Method 2 estimate')
  abline(0,1,lty=2,col='grey')
  box(lwd=1)
  points(m,
         m2,
         pch=21, 
         bg='grey')
  # arrows(m[m2>xylim[2]], 0.95*xylim[2], m[m2>xylim[2]], 0.99*xylim[2], length=0.02)
  # text(m[m2>xylim[2]], 0.93*xylim[2],round(m2[m2>xylim[2]],1),cex=0.7)
dev.off()     

