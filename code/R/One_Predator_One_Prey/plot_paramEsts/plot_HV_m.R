source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~

fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Hassell.Varley", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

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
  color.factor='Parasitoids', # 'None', 'Parasitoids' or 'Replacement'
  # color.vector=col.vec, # delete or specify above plot()
  xlab="Hassel-Varley interference rate (m)",
  ylab="",
  labels=TRUE,
  xlim=c(0,5)
)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
SS <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Hassell.Varley']]["50%",'exponent',"estimate"]))

parm <- exp(parm)

names(parm)<-sub('./Dataset_Code/','',names(parm))
names(parm)<-sub('.{2}$','',names(parm))
names(parm) <- paste0(names(parm), ' (',SS,')')

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_HV_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
ylim <- c(0,3)
plot(parm~SS,
     ylim=ylim,
     type='n',
     log='x', 
     xlab='Sample size (n)', 
     ylab='Hassel-Varley interference rate (m)')
abline(h=c(0,1),lty=2,col='grey')
arrows(SS[parm>ylim[2]], 1*ylim[2], SS[parm>ylim[2]], 1.03*ylim[2], length=0.02)
text(SS[parm>ylim[2]],0.98*ylim[2],round(parm[parm>ylim[2]],0),cex=0.5)
points(SS,parm,
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