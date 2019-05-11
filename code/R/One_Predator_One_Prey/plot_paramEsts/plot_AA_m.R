source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~

fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]


# data sets in which depletion occurred in over x% of replicates (non-replacement datasets)
depleted <- unlist(lapply(ffr.fits, depletion.check, cutoff=0.7))
# data sets on parasitoids
parasitoids <- unlist(lapply(ffr.fits, function(x){!x$study.info$predator}))

pch.vec<-rep(19,length(ffr.fits))
pch.vec[parasitoids] <- 1
pch.vec[depleted>0 & parasitoids] <- 0
pch.vec[depleted>0 & !parasitoids] <- 15

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')


pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits,
     model="Arditi.Akcakaya",
     parameter="exponent",
     ilink=exp,
     point.est='median',
     plot.SEs=FALSE,
     display.outlier.ests=TRUE,
     # color.factor='None', # 'None', 'Parasitoids' or 'Replacement'
     # color.vector=col.vec, # delete or specify above plot()
     # pch.factor='None', # 'None', 'Parasitoids' or 'Replacement'
     pch.vector=pch.vec, # delete or specify above plot()
     xlab="Arditi-Akcakaya interference rate (m)",
     ylab="",
     labels=labels,
     xlim=c(0,5)
  )
  legend('topright', 
         legend=c('Predator', 'Parasitoid', 'Replacement', 'Non-replacement', 'Profiled', 'Approximated','Bootstrapped'),
         pch=c(19,1,19,15,NA,NA,NA), 
         lty=c(NA,NA,NA,NA,'solid','dashed','dotted'),
         inset=0.05,
         cex=0.8)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))

parm <- exp(parm)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
ylim <- c(0,5)
plot(parm~sample.sizes,
     ylim=ylim,
     type='n',
     log='x', 
     xlab='Sample size (n)', 
     ylab='Arditi-Akcakaya interference rate (m)')
abline(h=c(0,1),lty=2,col='grey')
arrows(sample.sizes[parm>ylim[2]], 1*ylim[2], 
       sample.sizes[parm>ylim[2]], 1.03*ylim[2], 
       length=0.02)
text(sample.sizes[parm>ylim[2]],
     0.98*ylim[2],
     round(parm[parm>ylim[2]],0),
     cex=0.5)
points(sample.sizes, parm,
       pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm[parm<5],breaks=50,
       xlab='Arditi-Akcakaya interference rate (m)',
       main='',
       col='grey',
       ylim=c(0,16),
       xlim=c(0,5.5))
  abline(v=1,lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()


