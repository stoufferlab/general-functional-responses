source('../../lib/plot_coefs.R') # for order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

library(RColorBrewer)
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
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

# Grab summary of RMSD estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

RMSD.H1 <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Holling.I'][[1]][stat]}))
RMSD.H2 <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Holling.II'][[1]][stat]}))
RMSD.BD <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Beddington.DeAngelis'][[1]][stat]}))
RMSD.CM <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Crowley.Martin'][[1]][stat]}))
RMSD.SN1 <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Stouffer.Novak.I'][[1]][stat]}))

RMSDs <- data.frame(RMSD.H1, RMSD.H2, RMSD.BD, RMSD.CM, RMSD.SN1)
colnames(RMSDs) <- sub('RMSD.', '', colnames(RMSDs))
colnames(RMSDs)[5] <- "G"

# Define color of each model
CR<-brewer.pal(n = 8, name = 'RdBu')
Mcols <- c(CR[5:8],CR[4:1])
Mpch <- c(rep(21,5),rep(22,4))
  
minRMSDs <- apply(RMSDs, 1, min)
dRMSDs <- RMSDs - minRMSDs
rnkRMSDs <- t(apply(dRMSDs, 1, rank, ties.method='first'))
colnames(rnkRMSDs) <- colnames(RMSDs)

# Define delta RMSD cut-off for "well performing" models:
# Which models have an RMSD that is less than 1% of the data average (repeat for raw and bootstrapped datasets and merge)
delRMSDcutoff <- unlist(lapply(ffr.fits, function(x) 0.01*mean(x$study.info$data.Nconsumed)))
delRMSDcutoff2 <- unlist(lapply(ffr.fits, function(x) 0.01*mean(x$study.info$data.Nconsumed.mean)))
delRMSDcutoff[which(is.na(delRMSDcutoff))] <- delRMSDcutoff2[which(!is.na(delRMSDcutoff2))]

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_RMSD_ranks_dbs.pdf',height=6,width=2.25)
  par(mar=c(3,7,2.5,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
    plot(1:nrow(rnkRMSDs), 1:nrow(rnkRMSDs),
         type='n', yaxt='n',
         xlim=c(0.5,ncol(rnkRMSDs)+0.5),
         ylim=c(0,nrow(rnkRMSDs)+1),
         xlab='Model rank by RMSD',
         ylab='',
         axes=F)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white")
    axis(2, at=1:nrow(rnkRMSDs), labels=labels, cex.axis=0.5, las=2)
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0))

    # Which models have "reasonable" delta-RMSD?
    xats <-table(which(dRMSDs < delRMSDcutoff, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))
    # segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    # segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
    
    # shade behind ties
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.4
    polygon(pxats,pyats,col=grey(0.666),border=NA)
    
    for(m in 1:ncol(rnkRMSDs)){
      points(rnkRMSDs[,m], 1:nrow(rnkRMSDs), 
             type='p', col='black', 
             bg=Mcols[m], pch=Mpch[m],
             cex=1, lwd=0.2)
    }
    box(lwd=1)
    par(xpd=TRUE)
    legend(0.425,nrow(rnkRMSDs)+6,legend=colnames(rnkRMSDs),
           pch=Mpch,pt.bg=Mcols, col='black', bg='white',
           horiz=FALSE, pt.cex=1,cex=0.6, ncol=5, title='Model', bty='n')
dev.off()


# ~~~~~~~~~
# Count times each model is in each rank
Cnt<-apply(rnkRMSDs,2,function(x){table(factor(x,levels=1:ncol(rnkRMSDs)))})
Cnt
pCnt <- round(prop.table(Cnt,2)*100,1)
pCnt
# Concl:  AA is ranked 1st most frequently, followed by CM and BD.

# ~~~~~~~~~
# Either pass counts or combine counts and proportions before exporting table
tab_Cnt <- Cnt

# bCnt <- paste0(Cnt,' (',pCnt,')')
# bCnt[Cnt==0] <- 0
# tab_Cnt <- matrix(bCnt, nrow=nrow(Cnt), byrow=T, dimnames=dimnames(Cnt))

wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_figs/')

latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings.tex',label='table:RMSD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by RMSD.')

# latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings.tex',label='table:RMSD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) for which each functional response model achieved a given rank relative to all other models as judged by RMSD.')

setwd(wd)
stop()
# ~~~~~~~~~

# Of the times that AA is ranked first, how often are BD and CM  within X RMSD units?
cnt<-apply(dRMSDs[rnkRMSDs[,ncol(rnkRMSDs)]==1,3:4] < delRMSDcutoff, 2, sum)
cnt
Cnt[1,'AA']
cnt/Cnt[1,'AA']
# Concl: BD and CM are within cutoff 60% of the time.

# Of times that AA is ranked first, how often are *either* BD and CM within X RMSD units?
cnt<-sum(apply(dRMSDs[rnkRMSDs[,ncol(rnkRMSDs)]==1,3:4] < delRMSDcutoff, 1, sum)>0)
cnt
cnt/Cnt[1,'AA']
# Concl: BD or CM are within cutoff ~76% of the time

# ~~~~~~~~~
# What about for datasets that have a sample size of at least X?
SScut <- 50
Cnt<-apply(rnkRMSDs[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnkRMSDs)))})
Cnt
pCnt <- round(prop.table(Cnt,2)*100,1)
pCnt

# Either pass counts or combine counts and proportions before exporting table
tab_Cnt <- Cnt

# bCnt <- paste0(Cnt,' (',pCnt,')')
# bCnt[Cnt==0] <- 0
# tab_Cnt <- matrix(bCnt, nrow=nrow(Cnt), byrow=T, dimnames=dimnames(Cnt))

wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_figs/')

latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings_top50.tex',label='table:RMSD_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by RMSD.')

# latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings_top50.tex',label='table:RMSD_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by RMSD.')

setwd(wd)
# ~~~~~~~~~

