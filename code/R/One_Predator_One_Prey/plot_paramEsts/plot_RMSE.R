# Root mean square deviation (error)
# Error technically refers to out-of-sample
# Deviation technically refers to within-sample, hence use RMSD in manuscript

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

# Grab summary of RMSE estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

RMSE.H1 <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Holling.I'][[1]][stat]}))
RMSE.H2 <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Holling.II'][[1]][stat]}))
RMSE.BD <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Beddington.DeAngelis'][[1]][stat]}))
RMSE.CM <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Crowley.Martin'][[1]][stat]}))
RMSE.R <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Ratio'][[1]][stat]}))
RMSE.AG <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Arditi.Ginzburg'][[1]][stat]}))
RMSE.HV <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Hassell.Varley'][[1]][stat]}))
RMSE.AA <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Arditi.Akcakaya'][[1]][stat]}))

RMSEs <- data.frame(RMSE.H1, RMSE.H2, RMSE.BD, RMSE.CM, RMSE.R, RMSE.AG, RMSE.HV, RMSE.AA)
colnames(RMSEs) <- sub('RMSE.', '', colnames(RMSEs))

# Define color of each model
parm.k<-c(1,2,3,3, 1,2,2,3)
# CR<-brewer.pal(n = 8, name = 'RdBu')
# Mcols <- c(CR[5:8],CR[4:1])
# Mpch <- c(rep(21,4),rep(22,4))
CR<-brewer.pal(3,'YlOrRd')
Mcols <- c(CR[parm.k[1:4]],CR[parm.k[5:8]])
Mpch <- c(rep(21,4),rep(22,4))
Mpch2 <- c(NA,NA,NA,3,NA,NA,4,NA) # to overlay x's to differentiate models
  
minRMSEs <- apply(RMSEs, 1, min)
dRMSEs <- RMSEs - minRMSEs
rnkRMSEs <- t(apply(dRMSEs, 1, rank, ties.method='first'))
colnames(rnkRMSEs) <- colnames(RMSEs)

# Define delta RMSE cut-off for "well performing" models:
# Which models have an RMSE that is less than 1% of the data average (repeat for raw and bootstrapped datasets, which will each throw errors where the other does not apply, and merge)
delRMSEcutoff <- unlist(lapply(ffr.fits, 
                               function(x){0.01*mean(x$study.info$data.Nconsumed)}))
delRMSEcutoff2 <- unlist(lapply(ffr.fits, 
                                function(x){0.01*mean(x$study.info$data.Nconsumed.mean)}))
delRMSEcutoff[which(is.na(delRMSEcutoff))] <- delRMSEcutoff2[which(!is.na(delRMSEcutoff2))]

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_RMSE_ranks.pdf',height=6,width=2.25)
  par(mar=c(3,7,2.5,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
    plot(1:nrow(rnkRMSEs), 1:nrow(rnkRMSEs),
         type='n', yaxt='n',
         xlim=c(1,ncol(rnkRMSEs)),
         ylim=c(0,nrow(rnkRMSEs)+1),
         xlab='Model rank by RMSE',
         ylab='',
         axes=F)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white")
    axis(2, at=1:nrow(rnkRMSEs), labels=labels, cex.axis=0.5, las=2)
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0))

    # Which models have "reasonable" delta-RMSE?
    xats <-table(which(dRMSEs < delRMSEcutoff, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))+0.5
    segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
    
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.5
    polygon(pxats,pyats,col='grey90',border=NA)
    
    for(m in 1:ncol(rnkRMSEs)){
      points(rnkRMSEs[,m], 1:nrow(rnkRMSEs), 
             type='p', col='black', 
             bg=Mcols[m], pch=Mpch[m],
             cex=1, lwd=0.2)
      # Overlay x's
      points(rnkRMSEs[,m], 1:nrow(rnkRMSEs), 
             type='p',  col='black', 
             pch=Mpch2[m],
             cex=1, lwd=0.2)
    }
    box(lwd=1)
  par(xpd=TRUE)
    legend(-4,nrow(rnkRMSEs)+6,legend=colnames(rnkRMSEs),
           pch=Mpch,pt.bg=Mcols, col='black', bg='white',
           horiz=TRUE, pt.cex=1.1,cex=0.6, ncol=2, title='Model')
    # overlay x's
    legend(-4,nrow(rnkRMSEs)+6,legend=colnames(rnkRMSEs),
           pch=Mpch2, pt.bg=NA, col='black', bg=NA,
           horiz=TRUE, pt.cex=1.1,cex=0.6, ncol=2, title='')
dev.off()


# ~~~~~~~~~
# Count times each model is in each rank
Cnt<-apply(rnkRMSEs,2,function(x){table(factor(x,levels=1:ncol(rnkRMSEs)))})
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
setwd('../../../../results/R/OnePredOnePrey_tables/')

latex(tab_Cnt,file='OnePredOnePrey_RMSE_rankings.tex',label='table:RMSE_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by RMSE.')

# latex(tab_Cnt,file='OnePredOnePrey_RMSE_rankings.tex',label='table:RMSE_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) for which each functional response model achieved a given rank relative to all other models as judged by RMSE.')

setwd(wd)
# ~~~~~~~~~

# Of the times that AA is ranked first, how often are BD and CM  within X RMSE units?
cnt<-apply(dRMSEs[rnkRMSEs[,ncol(rnkRMSEs)]==1,3:4] < delRMSEcutoff, 2, sum)
cnt
Cnt[1,'AA']
cnt/Cnt[1,'AA']
# Concl: BD and CM are within cutoff 60% of the time.

# Of times that AA is ranked first, how often are *either* BD and CM within X RMSE units?
cnt<-sum(apply(dRMSEs[rnkRMSEs[,ncol(rnkRMSEs)]==1,3:4] < delRMSEcutoff, 1, sum)>0)
cnt
cnt/Cnt[1,'AA']
# Concl: BD or CM are within cutoff ~76% of the time

# ~~~~~~~~~
# What about for datasets that have a sample size of at least X?
SScut <- 50
Cnt<-apply(rnkRMSEs[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnkRMSEs)))})
Cnt
pCnt <- round(prop.table(Cnt,2)*100,1)
pCnt

# Either pass counts or combine counts and proportions before exporting table
tab_Cnt <- Cnt

# bCnt <- paste0(Cnt,' (',pCnt,')')
# bCnt[Cnt==0] <- 0
# tab_Cnt <- matrix(bCnt, nrow=nrow(Cnt), byrow=T, dimnames=dimnames(Cnt))

wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')

latex(tab_Cnt,file='OnePredOnePrey_RMSE_rankings_top50.tex',label='table:RMSE_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by RMSE.')

# latex(tab_Cnt,file='OnePredOnePrey_RMSE_rankings_top50.tex',label='table:RMSE_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by RMSE.')

setwd(wd)
# ~~~~~~~~~

