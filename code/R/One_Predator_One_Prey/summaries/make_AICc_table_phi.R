
# for the bundle_fits function
source('../../lib/plot_coefs.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grab summary of AICc estimates across bootstrapped fits
AICc.H1 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.I'][[1]][stat]}))
AICc.H2 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.II'][[1]][stat]}))
AICc.BD <- unlist(lapply(ffr.fits, function(x){ x$AICc['Beddington.DeAngelis'][[1]][stat]}))
AICc.CM <- unlist(lapply(ffr.fits, function(x){ x$AICc['Crowley.Martin'][[1]][stat]}))
AICc.SN1 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Stouffer.Novak.I'][[1]][stat]}))

# relabel the models
AICcs <- data.frame(AICc.H1, AICc.H2, AICc.BD, AICc.CM, AICc.SN1)
colnames(AICcs) <- sub('AICc.', '', colnames(AICcs))
colnames(AICcs)[5] <- "G"

# make things relative
minAICcs <- apply(AICcs, 1, min)
dAICcs <- AICcs - minAICcs

# get model ranks
rnkAICcs <- t(apply(dAICcs, 1, rank, ties.method='first'))
colnames(rnkAICcs) <- colnames(AICcs)
