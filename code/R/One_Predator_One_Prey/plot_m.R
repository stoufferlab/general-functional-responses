source('../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../lib/holling_method_one_predator_one_prey.R')
source('../lib/ratio_method_one_predator_one_prey.R')

load('../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# subset to include only studies considered by DeLong and Vasseur
# ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]

# fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="exponent")
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")

pdf('../../../results/R/OnePredOnePrey_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits[fit.order],
     model="Arditi.Akcakaya",
     parameter="exponent",
     ilink=exp,
     plot.SEs=TRUE,
     display.outlier.ests=TRUE,
     xlab="Arditi-Akcakaya interference rate (m)",
     ylab="",
     labels=TRUE,
     xlim=c(0,5)
  )
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
SS<-unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
m<-unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))
m<-exp(m)

names(m)<-sub('./Dataset_Code/','',names(m))
names(m)<-sub('.R.exponent','',names(m))
names(m) <- paste0(names(m), ' (',SS,')')


pdf('../../../results/R/OnePredOnePrey_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
plot(m~SS,
     ylim=c(0,5),
     type='n',
     log='x', 
     xlab='Sample size (n)', 
     ylab='Arditi-Akcakaya interference rate (m)')
abline(h=c(0,1),lty=2,col='grey')
points(SS,exp(m),
       pch=21, bg='grey')
dev.off()


pdf('../../../results/R/OnePredOnePrey_m_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(m[m<5],breaks=50,
       xlab='Arditi-Akcakaya interference rate (m)',
       main='',
       col='grey',
       ylim=c(0,16),
       xlim=c(0,5.5))
  abline(v=1,lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()


# To do:
# Plot of model AIC (ranks) by dataset
# Plot of DeLong estimates vs. sample size
# Plot of AAmethod2 estimates vs. sample size

