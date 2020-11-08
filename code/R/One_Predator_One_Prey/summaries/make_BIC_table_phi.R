
# for the bundle_fits function
source('../../lib/plot_coefs.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grab summary of BIC estimates across bootstrapped fits
BIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$BIC['Holling.I'][[1]][stat]}))
BIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$BIC['Holling.II'][[1]][stat]}))
BIC.BD <- unlist(lapply(ffr.fits, function(x){ x$BIC['Beddington.DeAngelis'][[1]][stat]}))
BIC.CM <- unlist(lapply(ffr.fits, function(x){ x$BIC['Crowley.Martin'][[1]][stat]}))
BIC.SN1 <- unlist(lapply(ffr.fits, function(x){ x$BIC['Stouffer.Novak.I'][[1]][stat]}))

# relabel the models
BICs <- data.frame(BIC.H1, BIC.H2, BIC.BD, BIC.CM, BIC.SN1)
colnames(BICs) <- sub('BIC.', '', colnames(BICs))
colnames(BICs)[5] <- "G"

# make things relative
minBICs <- apply(BICs, 1, min)
dBICs <- BICs - minBICs

# get model ranks
rnkBICs <- t(apply(dBICs, 1, rank, ties.method='first'))
colnames(rnkBICs) <- colnames(BICs)
