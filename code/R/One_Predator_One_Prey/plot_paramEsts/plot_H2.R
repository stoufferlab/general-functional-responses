source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ffr.fits <- profile_coefs(ffr.fits, 
                          model='Holling.II',
                          point.est='median',
                          printWarnings = TRUE,
                          which.pars = "attack")

save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.H2.Rdata')
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.H2.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General data and plot preparations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Holling.II", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

###################################################
# ~~~~~~~~~~~~~~~~~~ H2 attack ~ ~~~~~~~~~~~~~~~~~~
###################################################
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_H2_a.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
  model="Holling.II",
  parameter="attack",
  ilink=exp,
  point.est='median',
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
  xlab="Holling Type II attack rate (a)",
  labels=labels,
  xlim=c(0,1)
)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Alternative / additional summary plots
parm <- unlist(lapply(ffr.fits, function(x) x$estimates[['Holling.II']]["50%",'attack',"estimate"]))
parm <- exp(parm)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_H2_a_xy.pdf',height=3,width=4)
par(cex=0.7,  mgp=c(1.5,0.1,0), tcl=-0.1)
ylim <- c(1E-6,1E3)
plot(parm~sample.sizes,
     ylim=ylim,
     type='n',
     log='xy',
     xlab='Sample size (n)', 
     ylab='Holling Type II attack rate (a)')
arrows(sample.sizes[parm>ylim[2]], 1*ylim[2], sample.sizes[parm>ylim[2]], 1.03*ylim[2], length=0.02)
text(sample.sizes[parm>ylim[2]],0.98*ylim[2],round(parm[parm>ylim[2]],0), cex=0.5)
points(sample.sizes,parm,
       pch=21, bg='grey')
dev.off()

