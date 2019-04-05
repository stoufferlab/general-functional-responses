# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots of the raw data
# Future enhancement:  
# Grab data and fits from 'meta_analyze.R' output instead 
#      in order to add the fitted functions to the plots.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
  Sys.getenv("LOGNAME"),
  stouffer = '../../../dropbox_data/Data',
  marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# master list of datasets
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)

# remove template files which don't actually read data
datasets <- grep("template",datasets,invert=TRUE,value=TRUE)

# remove zzz files which are placeholders while a dataset is being cleaned/incorporated
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf('../../../results/R/OnePredOnePrey_DatasetPlots.pdf',height=10,width=6)
par(mfrow=c(8,4), mar=c(3,3,1.5,0.5), cex=0.6, cex.axis=0.8, cex.main=0.8, las=1, mgp=c(1.5,0.1,0), tcl=-0.1, pch=21)
for(i in 1:length(datasets)){
  message(datasets[i])
  source(datasets[i]) # loads the data into data frame 'd'
  
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
  
  title <- sub('./Dataset_Code/','',datasets[i])
  title <- sub('.R','',title)
  
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