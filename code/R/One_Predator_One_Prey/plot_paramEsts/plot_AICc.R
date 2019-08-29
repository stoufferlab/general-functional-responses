source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
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

# Grab summary of AICc estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

AICc.H1 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.I'][[1]][stat]}))
AICc.H2 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.II'][[1]][stat]}))
AICc.BD <- unlist(lapply(ffr.fits, function(x){ x$AICc['Beddington.DeAngelis'][[1]][stat]}))
AICc.CM <- unlist(lapply(ffr.fits, function(x){ x$AICc['Crowley.Martin'][[1]][stat]}))
AICc.R <- unlist(lapply(ffr.fits, function(x){ x$AICc['Ratio'][[1]][stat]}))
AICc.AG <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Ginzburg'][[1]][stat]}))
AICc.HV <- unlist(lapply(ffr.fits, function(x){ x$AICc['Hassell.Varley'][[1]][stat]}))
AICc.AA <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Akcakaya'][[1]][stat]}))

AICcs <- data.frame(AICc.H1, AICc.H2, AICc.BD, AICc.CM, AICc.R, AICc.AG, AICc.HV, AICc.AA)
colnames(AICcs) <- sub('AICc.', '', colnames(AICcs))

# Define color of each model
colnames(AICcs)
CR<-brewer.pal(n = 8, name = 'RdBu')
Mcols <- c(CR[5:8],CR[4:1])
Mpch <- c(rep(21,4),rep(22,4))

minAICcs <- apply(AICcs, 1, min)
dAICcs <- AICcs - minAICcs
rnkAICcs <- t(apply(dAICcs, 1, rank, ties.method='first'))
colnames(rnkAICcs) <- colnames(AICcs)

# Define delta AICc cut-off for "indistinguishably well performing" models
delAICcutoff <- 2

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_ranks.pdf',height=6,width=2.25)
  par(mar=c(3,7,2.5,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
    plot(1:nrow(rnkAICcs), 1:nrow(rnkAICcs),
         type='n', yaxt='n',
         xlim=c(1,ncol(rnkAICcs)),
         ylim=c(0,nrow(rnkAICcs)+1),
         xlab='Model rank by AICc',
         ylab='',
         axes=F)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
    axis(2, at=1:nrow(rnkAICcs), labels=labels, cex.axis=0.5, las=2)
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0))

    # Which models have delta-AICc within X=2 of best-performing model?
    xats <-table(which(dAICcs < delAICcutoff, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))+0.5
    segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
  
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.5
    polygon(pxats,pyats,col='grey90',border=NA)
    
    for(m in 1:ncol(rnkAICcs)){
      points(rnkAICcs[,m], 1:nrow(rnkAICcs), 
             type='p',  col='black', 
             bg=Mcols[m], pch=Mpch[m],
             cex=1, lwd=0.2)
    }  
    box(lwd=1)
  par(xpd=TRUE)
    legend(-4,nrow(rnkAICcs)+6,legend=colnames(rnkAICcs),
           pch=Mpch, pt.bg=Mcols, col='black', bg='white',
            horiz=TRUE, pt.cex=1.1,cex=0.6, ncol=2, title='Model')
dev.off()

# ~~~~~~~~~
# Count times each model is in each rank
Cnt<-apply(rnkAICcs,2,function(x){table(factor(x,levels=1:ncol(rnkAICcs)))})
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

latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings.tex',label='table:AICc_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

# latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings.tex',label='table:AICc_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

setwd(wd)
# ~~~~~~~~~

# Of the times that AA is ranked first, how often are BD and CM  within X AICc units?
cnt<-apply(dAICcs[rnkAICcs[,ncol(rnkAICcs)]==1,3:4] < delAICcutoff, 2, sum)
cnt
Cnt[1,'AA']
cnt/Cnt[1,'AA']
# Concl: BD is within 2 dAICc units > 1/2 the time, and CM is > 1/4 of the time

# Of times that AA is ranked first, how often are *either* BD and CM within X AICc units?
cnt<-sum(apply(dAICcs[rnkAICcs[,ncol(rnkAICcs)]==1,3:4] < delAICcutoff, 1, sum)>0)
cnt
cnt/Cnt[1,'AA']
# Concl: BD or CM are within 2 dAICc units ~60% of the time

# ~~~~~~~~~
# What about for datasets that have a sample size of at least X?
SScut <- 50
Cnt<-apply(rnkAICcs[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnkAICcs)))})
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

latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings_top50.tex',label='table:AICc_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

# latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings_top50.tex',label='table:AICc_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency)  with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

setwd(wd)
# ~~~~~~~~~

#~~~~~~
# dAICc
#~~~~~~
# plot(1:nrow(dAICcs), 1:nrow(dAICcs),
#      type='n', yaxt='n',
#      # xlim=c(0,max(dAICcs)),
#      xlim=c(0,4),
#      xlab='Model delta-AICc',
#      ylab=''
# )
# axis(side=2, at=1:nrow(dAICcs), labels=labels, cex.axis=0.5, las=2)
# 
# Mcols <- rainbow(ncol(dAICcs))
# for(m in 1:ncol(dAICcs)){
#   points(dAICcs[,m], 1:nrow(dAICcs), pch=21, col=Mcols[m], bg=Mcols[m])
# }  
# 
# legend('bottomright',legend=colnames(dAICcs),pch=21, pt.bg=Mcols, col=Mcols)
# 
