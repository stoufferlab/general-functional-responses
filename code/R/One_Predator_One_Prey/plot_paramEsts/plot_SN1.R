source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')
# source('../../lib/Bhat/dfp.R')
# source('../../lib/Bhat/dqstep.R')
# source('../../lib/Bhat/logit.hessian.R')
# source('../../lib/Bhat/plkhci.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# to color models that are not supported by AICc differently
stat <- '50%'
# only include fits for which the SN1 model is at least tied based on AICc
AICc.H1 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.I'][[1]][stat]}))
AICc.H2 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.II'][[1]][stat]}))
AICc.BD <- unlist(lapply(ffr.fits, function(x){ x$AICc['Beddington.DeAngelis'][[1]][stat]}))
AICc.CM <- unlist(lapply(ffr.fits, function(x){ x$AICc['Crowley.Martin'][[1]][stat]}))
AICc.SN1 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Stouffer.Novak.I'][[1]][stat]}))
# AICc.R <- unlist(lapply(ffr.fits, function(x){ x$AICc['Ratio'][[1]][stat]}))
# AICc.AG <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Ginzburg'][[1]][stat]}))
# AICc.HV <- unlist(lapply(ffr.fits, function(x){ x$AICc['Hassell.Varley'][[1]][stat]}))
# AICc.AA <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Akcakaya'][[1]][stat]}))

AICcs <- data.frame(AICc.H1, AICc.H2, AICc.BD, AICc.CM, AICc.SN1)
colnames(AICcs) <- sub('AICc.', '', colnames(AICcs))
colnames(AICcs)[5] <- "G"
minAICcs <- apply(AICcs, 1, min)
dAICcs <- AICcs - minAICcs

# only include "supported" model fits
# ffr.fits <- ffr.fits[which(dAICcs$G <= 2)]

color.vector <- numeric(length(ffr.fits))
for(i in 1:length(ffr.fits)){
  if(dAICcs$G[i] > 2){
    color.vector[i] <- "red"
  }else{
    color.vector[i] <- "black"
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ffr.fits <- profile_coefs(ffr.fits, 
#                           model='Stouffer.Novak.I',
#                           point.est='median',
#                           printWarnings = TRUE)

# save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')
load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General data and plot preparations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Stouffer.Novak.I", order.parm="phi_denom")
ffr.fits <- ffr.fits[fit.order]

color.vector <- color.vector[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

###################################################
# ~~~~~~~~~~~~~~~~~~ SN1 PhiDenom ~~~~~~~~~~~~~~~~~
###################################################
pdf(file="../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_SN1_PhiDenom.pdf",height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7)
plot.coefs(
  ffr.fits,
	model="Stouffer.Novak.I",
	parameter="phi_denom",
  ilink=identity,
  point.est='median',
  plot.SEs=TRUE,
  display.outlier.ests=TRUE,
	xlab="Effect of feeding on interfering",
  labels=labels,
  vertLines = c(0,1),
	xlim=c(-3,3),
  color.vector = color.vector
)
dev.off()

