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

# Grab summary of RMSD estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

RMSD.H1 <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Holling.I'][[1]][stat]}))
RMSD.H2 <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Holling.II'][[1]][stat]}))
RMSD.BD <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Beddington.DeAngelis'][[1]][stat]}))
RMSD.CM <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Crowley.Martin'][[1]][stat]}))
RMSD.R <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Ratio'][[1]][stat]}))
RMSD.AG <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Arditi.Ginzburg'][[1]][stat]}))
RMSD.HV <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Hassell.Varley'][[1]][stat]}))
RMSD.AA <- unlist(lapply(ffr.fits, function(x){ x$RMSE['Arditi.Akcakaya'][[1]][stat]}))

RMSDs <- data.frame(RMSD.H1, RMSD.H2, RMSD.BD, RMSD.CM, RMSD.R, RMSD.AG, RMSD.HV, RMSD.AA)
colnames(RMSDs) <- sub('RMSD.', '', colnames(RMSDs))

# Define color of each model
colnames(RMSDs)
parm.k.1<-c(1,2,3,3)
parm.k.2<-c(1,2,2,3)
CR1<-brewer.pal(n = 3, name = 'YlOrRd')
CR2<-brewer.pal(n = 3, name = 'Blues')
Mcols<-c(CR1[parm.k.1],CR2[parm.k.2])

Mpch <- c(rep(21,4),rep(22,4))
Mpch2 <- c(NA,NA,NA,21,NA,NA,22,NA) # overlay symbols to differentiate equal k models

minRMSDs <- apply(RMSDs, 1, min)
delV.RMSD <- RMSDs - minRMSDs
rnks.RMSD <- t(apply(delV.RMSD, 1, rank, ties.method='first'))
colnames(rnks.RMSD) <- colnames(RMSDs)

# Define delta RMSD cut-off for "well performing" models:
# Which models have an RMSD that is less than 1% of the data average (repeat for raw and bootstrapped datasets, which will each throw errors where the other does not apply, and merge)
cut.RMSD <- unlist(lapply(ffr.fits, 
                               function(x){0.01*mean(x$study.info$data.Nconsumed)}))
cut.RMSD2 <- unlist(lapply(ffr.fits, 
                                function(x){0.01*mean(x$study.info$data.Nconsumed.mean)}))
cut.RMSD[which(is.na(cut.RMSD))] <- cut.RMSD2[which(!is.na(cut.RMSD2))]

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_RMSD_ranks.pdf',
    height=6.5,width=2)
par(mar=c(2,6.5,4,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
    plot(1:nrow(rnks.RMSD), 1:nrow(rnks.RMSD),
         type='n', yaxt='n',
         xlim=c(1,ncol(rnks.RMSD)),
         ylim=c(0,nrow(rnks.RMSD)+1),
         xlab='',
         ylab='',
         axes=F)
    title(xlab='Model rank by RMSD',line=1)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white")
    axis(2, at=1:nrow(rnks.RMSD), labels=labels, cex.axis=0.5, las=2)
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0))

    # Which models have "reasonable" delta-RMSD?
    xats <-table(which(delV.RMSD < cut.RMSD, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))+0.5
    segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
    
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.5
    polygon(pxats,pyats,col='grey90',border=NA)
    
    for(m in 1:ncol(rnks.RMSD)){
      points(rnks.RMSD[,m], 1:nrow(rnks.RMSD), 
             type='p',  col='black', 
             bg=Mcols[m], pch=Mpch[m],
             cex=1, lwd=0.2)
      # Overlay secondary symbos
      points(rnks.RMSD[,m], 1:nrow(rnks.RMSD), 
             type='p',  col='black', 
             bg='white',pch=Mpch2[m],
             cex=0.3, lwd=0.2)
    } 
    box(lwd=1)
    
    par(xpd=TRUE)
    vadj <- 7.5
    xadj <- 1.5
    mods<-c(1,5)
    legend(xadj,nrow(rnks.RMSD)+vadj,legend=colnames(rnks.RMSD)[mods],
           pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
           pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8, 
           title=expression(italic(k)==1), bty='n')
    # overlay secondary symbols
    legend(xadj,nrow(rnks.RMSD)+vadj,legend=colnames(rnks.RMSD)[mods],
           pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
           pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
           title=expression(italic(k)==1), bty='n')
    
    mods<-c(2,6,7)
    legend(xadj+2,nrow(rnks.RMSD)+vadj,legend=colnames(rnks.RMSD)[mods],
           pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
           pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
           title=expression(italic(k)==2), bty='n')
    # overlay secondary symbols
    legend(xadj+2,nrow(rnks.RMSD)+vadj,legend=colnames(rnks.RMSD)[mods],
           pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
           pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
           title=expression(italic(k)==2), bty='n')
    
    mods<-c(3,4,8)
    legend(xadj+4,nrow(rnks.RMSD)+vadj,legend=colnames(rnks.RMSD)[mods],
           pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
           pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
           title=expression(italic(k)==3), bty='n')
    # overlay secondary symbols
    legend(xadj+4,nrow(rnks.RMSD)+vadj,legend=colnames(rnks.RMSD)[mods],
           pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
           pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
           title=expression(italic(k)==3), bty='n')
    
    rect(1.5,nrow(rnks.RMSD)+vadj+1.75,7.5,nrow(rnks.RMSD)+vadj-6)
    text(xadj+3,nrow(rnks.RMSD)+vadj+0.75,'Models', adj=0.5, cex=0.8)
dev.off()


# ~~~~~~~~~
# Count times each model is in each rank
Cnt<-apply(rnks.RMSD,2,function(x){table(factor(x,levels=1:ncol(rnks.RMSD)))})
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

latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings.tex',label='table:RMSD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by the root mean square deviation (RMSD).')

# latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings.tex',label='table:RMSD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) for which each functional response model achieved a given rank relative to all other models as judged by RMSD.')

setwd(wd)
# ~~~~~~~~~

# Of the times that AA is ranked first, how often are BD and CM  within X RMSD units?
cnt<-apply(delV.RMSD[rnks.RMSD[,ncol(rnks.RMSD)]==1,3:4] < cut.RMSD, 2, sum)
cnt
Cnt[1,'AA']
cnt/Cnt[1,'AA']
# Concl: BD and CM are within cutoff 60% of the time.

# Of times that AA is ranked first, how often are *either* BD and CM within X RMSD units?
cnt<-sum(apply(delV.RMSD[rnks.RMSD[,ncol(rnks.RMSD)]==1,3:4] < cut.RMSD, 1, sum)>0)
cnt
cnt/Cnt[1,'AA']
# Concl: BD or CM are within cutoff ~76% of the time

# ~~~~~~~~~
# What about for datasets that have a sample size of at least X?
SScuts <- seq(5,300,by=1)
fFirst.RMSD<-fSecnd.RMSD<-dim(0)
for(SScut in SScuts){
  Cnt<-apply(rnks.RMSD[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnks.RMSD)))})
  pCnt <- prop.table(Cnt,2)
  fFirst.RMSD <- rbind(fFirst.RMSD,pCnt[1,])
  fSecnd.RMSD <- rbind(fSecnd.RMSD,pCnt[2,])
}

# Darken colors of linear models
Mcols[c(1,5)]<-'grey40'
ltys <- c(1,1,1,6,2,2,3,1)
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_RMSD_rankBySS.pdf',
    height=4,width=3)
par(mar=c(2,6.5,4,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.8)
layout(matrix(c(1,2,3),ncol=1), heights=c(1,3,3)) 
par(mar=c(0,3,0,0))
plot(1,1,type="n", axes=F, xlab="", ylab="")

xadj<-0.7
yadj<-1.2
lcex<-0.7

mods<-c(1,5)
legend(xadj,yadj,legend=colnames(fFirst.RMSD)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       cex=lcex, ncol=1, pt.lwd=0.8, 
       title=expression(italic(k)==1), bty='n')

mods<-c(2,6,7)
legend(xadj+0.2,yadj,legend=colnames(fFirst.RMSD)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==2), bty='n')

mods<-c(3,4,8)
legend(xadj+0.4,yadj,legend=colnames(fFirst.RMSD)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==3), bty='n')

text(xadj+0.3, yadj+0.1, 'Models', adj=0.5, cex=1)


par(mar=c(3,3,0.5,1))
matplot(SScuts,fFirst.RMSD,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
        ylim=c(0,0.8),lwd=1.5,
        ylab='Fraction in first rank',
        xlab='Sample size greater than...')
legend('topleft',legend='A',inset=-0.05,bty='n',cex=1.2)
matplot(SScuts,fSecnd.RMSD,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
        ylim=c(0,0.4),lwd=1.5,
        ylab='Fraction in second rank',
        xlab='Sample size greater than...')
legend('topleft',legend='B',inset=-0.05,bty='n',cex=1.2)
dev.off()



SScut <- 50
Cnt<-apply(rnks.RMSD[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnks.RMSD)))})
Cnt
sum(Cnt[,1])
pCnt <- prop.table(Cnt,2)
pCnt

# Either pass counts or combine counts and proportions before exporting table
tab_Cnt <- Cnt

# bCnt <- paste0(Cnt,' (',pCnt,')')
# bCnt[Cnt==0] <- 0
# tab_Cnt <- matrix(bCnt, nrow=nrow(Cnt), byrow=T, dimnames=dimnames(Cnt))

wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')

latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings_top50.tex',label='table:RMSD_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by the root mean square deviation (RMSD).')

# latex(tab_Cnt,file='OnePredOnePrey_RMSD_rankings_top50.tex',label='table:RMSD_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by RMSD.')

setwd(wd)
# ~~~~~~~~~

