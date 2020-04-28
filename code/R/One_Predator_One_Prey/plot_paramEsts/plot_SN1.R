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
AIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.I'][[1]][stat]}))
AIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.II'][[1]][stat]}))
AIC.BD <- unlist(lapply(ffr.fits, function(x){ x$AIC['Beddington.DeAngelis'][[1]][stat]}))
AIC.CM <- unlist(lapply(ffr.fits, function(x){ x$AIC['Crowley.Martin'][[1]][stat]}))
AIC.SN1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Stouffer.Novak.I'][[1]][stat]}))
# AICc.R <- unlist(lapply(ffr.fits, function(x){ x$AICc['Ratio'][[1]][stat]}))
# AICc.AG <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Ginzburg'][[1]][stat]}))
# AICc.HV <- unlist(lapply(ffr.fits, function(x){ x$AICc['Hassell.Varley'][[1]][stat]}))
# AICc.AA <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Akcakaya'][[1]][stat]}))

AICs <- data.frame(AIC.H1, AIC.H2, AIC.BD, AIC.CM, AIC.SN1)
colnames(AICs) <- sub('AIC.', '', colnames(AICs))
colnames(AICs)[5] <- "G"
minAICs <- apply(AICs, 1, min)
dAICs <- AICs - minAICs

# fits for which the SN1 model is not supported are colored red
color.vector <- numeric(length(ffr.fits))
for(i in 1:length(ffr.fits)){
  if(dAICs$G[i] > 2){
    color.vector[i] <- "red"
  }else{
    color.vector[i] <- "black"
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ffr.fits <- profile_coefs(
  ffr.fits,
  model='Stouffer.Novak.I',
  point.est='median',
  printWarnings = FALSE,
  which.pars = "phi_denom"
)
save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

