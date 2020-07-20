# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots of the raw data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_DatasetPlots.pdf',height=10,width=6)
par(mfrow=c(8,4), mar=c(3,3,1.5,0.5), cex=0.6, cex.axis=0.8, cex.main=0.8, las=1, mgp=c(1.5,0.1,0), tcl=-0.1, pch=21)
for(i in 1:length(ffr.fits)){
  message(ffr.fits[[i]]$study.info$datasetName)
  
  # figure out which columns pertain to data
  cols <- names(ffr.fits[[i]]$study.info)
  dcols <- cols[grepl("data[.]", cols)]
  d <- as.data.frame(ffr.fits[[i]]$study.info[dcols])
  colnames(d) <- gsub("data[.]","",colnames(d))
  
  #Create a function to generate a continuous color palette for predator densities
  rbPal <- colorRampPalette(c('red','blue'))
  d$Col <- rbPal(10)[as.numeric(cut(d$Npredator, breaks = 10))]
  
  d$Ratio <- d$Nprey/d$Npredator
  if("Nconsumed.mean" %in% colnames(d)){
    d$Nconsumed.pPred.mean <- d$Nconsumed.mean/d$Npredator
    d$Nconsumed.pPred.se <- d$Nconsumed.se/d$Npredator
  }else{
    d$Nconsumed.pPred <- d$Nconsumed/d$Npredator
  }
  
  title <- ffr.fits[[i]]$study.info$datasetName
  
  if("Nconsumed.mean" %in% colnames(d)){
    # ~~~~~~ Total eaten
    ylims <- c(0, max(d$Nconsumed.mean+d$Nconsumed.se))
    if(is.na(ylims[2])){
      ylims[2] <- max(d$Nconsumed.mean)
    }
    xlims <- c(0, max(d$Nprey))
    plot(d$Nconsumed.mean~d$Nprey, 
         xlab='Prey', ylab='Eaten',
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    suppressWarnings(
      arrows(d$Nprey,d$Nconsumed.mean-d$Nconsumed.se, 
             d$Nprey,d$Nconsumed.mean+d$Nconsumed.se, 
             length=0.05, angle=90, code=3)      )
    
    xlims <- c(0, max(d$Ratio))
    plot(d$Nconsumed.mean~d$Ratio, 
         xlab='Prey/Pred', ylab='Eaten', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    suppressWarnings(
      arrows(d$Ratio,d$Nconsumed.mean-d$Nconsumed.se, d$Ratio,d$Nconsumed.mean+d$Nconsumed.se, 
             length=0.05, angle=90, code=3)    )
    
    # ~~~~~~ Eaten per predator
    xlims <- c(0, max(d$Nprey))
    plot(d$Nconsumed.pPred.mean~d$Nprey, 
         xlab='Prey', ylab='Eaten/Pred', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    suppressWarnings(
      arrows(d$Nprey,d$Nconsumed.pPred.mean-d$Nconsumed.pPred.se, d$Nprey,d$Nconsumed.pPred.mean+d$Nconsumed.pPred.se,
             length=0.05, angle=90, code=3)      )
    
    xlims <- c(0, max(d$Ratio))
    plot(d$Nconsumed.pPred.mean~d$Ratio, 
         xlab='Prey/Pred', ylab='Eaten/Pred', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    suppressWarnings(
      arrows(d$Ratio,d$Nconsumed.pPred.mean-d$Nconsumed.pPred.se, d$Ratio,d$Nconsumed.pPred.mean+d$Nconsumed.pPred.se,
             length=0.05, angle=90, code=3)      )
  }else{
    # ~~~~~~ Total eaten
    ylims <- c(0, max(d$Nconsumed))
    
    xlims <- c(0, max(d$Nprey))
    plot(d$Nconsumed~d$Nprey, 
         xlab='Prey', ylab='Eaten', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    
    xlims <- c(0, max(d$Ratio))
    plot(d$Nconsumed~d$Ratio, 
         xlab='Prey/Pred', ylab='Eaten', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    
    # ~~~~~~ Eaten per predator
    ylims <- c(0, max(d$Nconsumed.pPred))
    
    xlims <- c(0, max(d$Nprey))
    plot(d$Nconsumed.pPred~d$Nprey, 
         xlab='Prey', ylab='Eaten/Pred', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
    
    xlims <- c(0, max(d$Ratio))
    plot(d$Nconsumed.pPred~d$Ratio, 
         xlab='Prey/Pred', ylab='Eaten/Pred', 
         bg=d$Col, xlim=xlims, ylim=ylims, main=title)
  }
}
dev.off()