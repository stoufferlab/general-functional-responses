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
ffr.fits <- profile_coefs(
  ffr.fits,
  model='Arditi.Akcakaya',
  point.est='median',
  printWarnings = FALSE,
  which.pars=c("attack","exponent")
)
save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.AA.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.AA.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General data and plot preparations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

# parasitoid data sets
parasitoids <- unlist(lapply(ffr.fits, function(x){!x$study.info$predator}))

# replacement studies
replacement <- unlist(lapply(ffr.fits, function(x){!x$study.info$replacement}))

# data sets in which depletion (90% of prey available prey consumed) occurred in over 10% of replicates
depleted <- unlist(lapply(ffr.fits, depletion.check))

pch.vec<-rep(19,length(ffr.fits))
pch.vec[parasitoids] <- 15
pch.vec[!parasitoids & !replacement] <- 1
pch.vec[!parasitoids & !replacement & depleted] <- 10
pch.vec[parasitoids & !replacement] <- 0
pch.vec[parasitoids & !replacement & depleted] <- 7

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

###################################################
# ~~~~~~~~~~~~~~~~~~ AA exponent ~~~~~~~~~~~~~~~~~~
###################################################
cairo_pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
  model="Arditi.Akcakaya",
  parameter="exponent",
  ilink=exp,
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
  xlab=expression(paste("Interference rate ",(italic(m)))),
  labels=labels,
  xlim=c(0,5),
  vertLines=c(0,1)
)
dev.off()

cairo_pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m_v2.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
  plot.coefs(
     ffr.fits,
     model="Arditi.Akcakaya",
     parameter="exponent",
     ilink=exp,
     plot.SEs=TRUE,
     SE.lty=c(1,2,3),
     display.outlier.ests=TRUE,
     pch.vector=pch.vec,
     xlab=expression(paste("Interference rate ",(italic(m)))),
     labels=labels,
     xlim=c(0,5),
     vertLines=c(0,1)
  )
  legend(3,70, 
         legend=c('Predator', 'Parasitoid', NA,
                  'Replacement', 'Non-replacement', NA,
                  'No depletion','Depletion',NA,
                  'Profiled', 'Approximated','Bootstrapped'),
         pch=c(19,15,NA,19,1,NA,1,10,NA,NA,NA,NA), 
         lty=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,1,2,3),
         inset=0.05)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))

out <- data.frame(dataset=labels,m=round(parm,3))
write.csv(out,
          '../../../../results/R/OnePredOnePrey_tables/OnePredOnePrey_AA-m_estimates.csv',
          row.names=FALSE)

parm <- exp(parm)

# Summary stats
round(range(parm),2)
round(median(parm),2)
round(mean(parm),2)
round(sd(parm),2)

# Summary stats - for studies exceeding median sample size
s <- median(sample.sizes)
round(range(parm[sample.sizes>s]),2)
round(median(parm[sample.sizes>s]),2)
round(mean(parm[sample.sizes>s]),2)
round(sd(parm[sample.sizes>s]),2)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  plot(parm~sample.sizes,
       type='n',
       log='x', 
       xlab='Sample size (n)', 
       ylab='Arditi-Akcakaya interference rate (m)')
  abline(h=c(0,1),lty=2,col='grey')
  points(sample.sizes, parm,
         pch=21, bg='grey')
dev.off()

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm[parm<5],breaks=50,
       xlab='Arditi-Akcakaya interference rate (m)',
       main='',
       col='grey')
  abline(v=1,lty=3,lwd=1,col='black')
  box(lwd=1,bty='l')
dev.off()

######################################################
# ~~~~~~~~~~~~~~~~~~ AA attack rate ~~~~~~~~~~~~~~~~~~
######################################################
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_a.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
  model="Arditi.Akcakaya",
  parameter="attack",
  ilink=exp,
  point.est='median',
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
  pch.vector=pch.vec,
  xlab=expression(paste("Attack rate ",(italic(a)))),
  labels=labels,
  xlim=c(0,1)
)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'attack',"estimate"]))

parm <- exp(parm)

# Summary stats
round(range(parm),2)
round(median(parm),2)
round(mean(parm),2)
round(sd(parm),2)

# Summary stats - for studies exceeding median sample size
s <- median(sample.sizes)
round(range(parm[sample.sizes>s]),2)
round(median(parm[sample.sizes>s]),2)
round(mean(parm[sample.sizes>s]),2)
round(sd(parm[sample.sizes>s]),2)

names(parm)<-sub('./Dataset_Code/','',names(parm))
names(parm)<-sub('.{2}$','',names(parm))
names(parm) <- paste0(names(parm), ' (',sample.sizes,')')

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_a_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  ylim <- c(0,5)
  plot(parm~sample.sizes,
       ylim=ylim,
       type='n',
       log='x', 
       xlab='Sample size (n)', 
       ylab='Arditi-Akcakaya attack rate (a)')
  arrows(sample.sizes[parm>ylim[2]], 1*ylim[2], sample.sizes[parm>ylim[2]], 1.03*ylim[2], length=0.02)
  text(sample.sizes[parm>ylim[2]],0.98*ylim[2],round(parm[parm>ylim[2]],0),cex=0.5)
  points(sample.sizes,parm,
         pch=21, bg='grey')
dev.off()


pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_a_hist.pdf',height=2.5,width=4)
par(mar=c(3,3,1,1), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i',xaxs='i')
  hist(parm[parm<5],breaks=50,
       xlab='Arditi-Akcakaya attack rate (a)',
       main='',
       col='grey',
       ylim=c(0,60),
       xlim=c(0,5.5))
  box(lwd=1,bty='l')
dev.off()

