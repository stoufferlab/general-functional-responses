source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

load('../../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~

fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

# fraction of replicates in which significant depletion occurred (in non-replacement datasets)
col.vec<-rep('black',length(ffr.fits))
depleted <- unlist(lapply(ffr.fits, depletion.check, cutoff=0.7))
col.vec[depleted>0] <- 'red'

pdf('../../../../results/R/OnePredOnePrey_AA_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits,
     model="Arditi.Akcakaya",
     parameter="exponent",
     ilink=exp,
     point.est='median',
     plot.SEs=TRUE,
     display.outlier.ests=TRUE,
     color.factor='None', # 'None', 'Parasitoids' or 'Replacement'
     # color.vector=col.vec, # delete or specify above plot()
     pch.factor='Parasitoids', # 'None', 'Parasitoids' or 'Replacement'
     # pch.vector=col.vec, # delete or specify above plot()
     xlab="Arditi-Akcakaya interference rate (m)",
     ylab="",
     labels=TRUE,
     xlim=c(0,5)
  )
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
SS <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))

parm <- exp(parm)

names(parm)<-sub('./Dataset_Code/','',names(parm))
names(parm)<-sub('.{2}$','',names(parm))
names(parm) <- paste0(names(parm), ' (',SS,')')

pdf('../../../../results/R/OnePredOnePrey_AA_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
ylim <- c(0,5)
plot(parm~SS,
     ylim=ylim,
     type='n',
     log='x', 
     xlab='Sample size (n)', 
     ylab='Arditi-Akcakaya interference rate (m)')
abline(h=c(0,1),lty=2,col='grey')
arrows(SS[parm>ylim[2]], 1*ylim[2], SS[parm>ylim[2]], 1.03*ylim[2], length=0.02)
text(SS[parm>ylim[2]],0.98*ylim[2],round(parm[parm>ylim[2]],0),cex=0.5)
points(SS,parm,
       pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_AA_m_hist.pdf',height=2.5,width=4)
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


