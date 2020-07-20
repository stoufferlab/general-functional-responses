source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Studies considered by DeLong and Vasseur
ffr.fits <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

m <- unlist(lapply(ffr.fits, function(x) x$estimates[['Arditi.Akcakaya']]["50%",'exponent',"estimate"]))

dat <- read.csv('Compilation_MLEs/DeLong_MLEs.csv')

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AA_m_xy_DeLong.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
  plot(dat$nObs,
       dat$DeLongM,
       type='n',
       xlab='Sample size (n)',
       ylab='DeLong & Vasseur\nArditi-Akcakaya interference rate (m)')
  abline(h=c(0,1),lty=2,col='grey')
  points(dat$nObs,
         dat$DeLongM,
         pch=21, bg='grey')
dev.off()
