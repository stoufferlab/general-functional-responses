source('../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../lib/plot_obsVfit.R')
source('../lib/holling_method_one_predator_one_prey.R')
source('../lib/ratio_method_one_predator_one_prey.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('../../../results/R/OnePredOnePrey_ffr.fits.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[rev(fit.order)]

n <- length(ffr.fits)
models <- c('Holling.I','Holling.II','Ratio','Arditi.Akcakaya','Beddington.DeAngelis')

pdf(file="../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_obsVfit.pdf",height=2.5,width=8, onefile = T)
par(mar=c(2,3,0.5,0.5), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.4)
par(mfcol=c(1, length(models)))
for(i in 1:length(ffr.fits)){
  titles <- models
  dataset <- ffr.fits[[i]]$study.info$datasetName
  titles[1] <- paste(dataset,'\n', titles[1])
  for(m in 1:length(models)){
    plot_obsVfit(ffr.fits[[i]], models[m], title=titles[m])
    mtext(bquote(LL==
                   .(round(ffr.fits[[i]]$LL[[models[m]]]['50%'],1)) 
                    ~ "("*.(round(ffr.fits[[i]]$LL[[models[m]]]['16%'],1))
                    *","~.(round(ffr.fits[[i]]$LL[[models[m]]]['84%'],1))*")"),
          side=1,line=3,cex=0.5)
    mtext(bquote(AIC==
                   .(round(ffr.fits[[i]]$AIC[[models[m]]]['50%'],1)) 
                 ~ "("*.(round(ffr.fits[[i]]$AIC[[models[m]]]['16%'],1))
                 *","~.(round(ffr.fits[[i]]$AIC[[models[m]]]['84%'],1))*")"),
          side=1,line=4,cex=0.5)
  }}
dev.off()
