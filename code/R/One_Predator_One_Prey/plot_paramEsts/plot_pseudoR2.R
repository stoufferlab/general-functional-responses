source('../../lib/plot_coefs.R') # for order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
load('../../../../results/R/OnePredOnePrey_ffr.fits.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

models <- c('Holling.Type.I','Holling.Type.II','Ratio','Arditi.Akcakaya','Beddington.DeAngelis')

