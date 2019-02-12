source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

load('../../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove datasets where AA2 method was NOT successfully applied
noAA2 <- unlist(lapply(ffr.fits, function(x) is.na(x$estimates['Arditi.Akcakaya.Method.2'])))
ffr.fits[which(noAA2)] <- NULL

fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya.Method.2", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

pdf('../../../../results/R/OnePredOnePrey_AA2_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits,
     model="Arditi.Akcakaya.Method.2",
     parameter="exponent",
     ilink=exp,
     plot.SEs=TRUE,
     point.est='median',
     display.outlier.ests=TRUE,
     xlab="Arditi-Akcakaya interference rate (m) (Method 2)",
     ylab="",
     labels=TRUE
  )
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Additional summary plots
SS <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
m2 <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya.Method.2']]["50%",'exponent',"estimate"]))

parm <- exp(m2)

names(parm)<-sub('./Dataset_Code/','',names(parm))
names(parm)<-sub('.R.exponent','',names(parm))
names(parm) <- paste0(names(parm), ' (',SS,')')

pdf('../../../../results/R/OnePredOnePrey_AA2_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
ylim <- c(0,1.7)
plot(parm~SS,
     ylim=ylim,
     type='n',
     log='x', 
     xlab='Sample size (n)', 
     ylab='Arditi-Akcakaya interference rate (m) (AA2)')
abline(h=c(0,1),lty=2,col='grey')
points(SS,parm,
       pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_AA2_m_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm[parm<5],breaks=50,
       xlab='Arditi-Akcakaya interference rate (m) (AA2)',
       main='',
       col='grey',
       ylim=c(0,6.5),
       xlim=c(0,1.7))
  abline(v=1,lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()

# Compare AA2 estimates to those of non-AA2 method
m <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))

pdf('../../../../results/R/OnePredOnePrey_AA2_m_vs_AA_m.pdf',height=3,width=3)
par(pty='s',mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7)
plot(exp(m),exp(m2),
     type='n',
     xlab='Maximum likelihood estimate',
     ylab='AA Method 2 estimate')
abline(0,1,lty=2,col='grey')
points(exp(m),
       exp(m2),
       pch=21, 
       bg='grey')
dev.off()     
