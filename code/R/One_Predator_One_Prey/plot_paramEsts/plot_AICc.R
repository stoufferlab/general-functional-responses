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
AICc.HV <- unlist(lapply(ffr.fits, function(x){ x$AICc['Hassell.Varley'][[1]][stat]}))
AICc.AG <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Ginzburg'][[1]][stat]}))
AICc.AA <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Akcakaya'][[1]][stat]}))

AICcs <- data.frame(AICc.H1, AICc.H2, AICc.BD, AICc.CM, AICc.R, AICc.HV, AICc.AG, AICc.AA)
colnames(AICcs) <- sub('AICc.', '', colnames(AICcs))

# Define color of each model
colnames(AICcs)
parm.k.1<-c(1,2,3,3)
parm.k.2<-c(1,2,2,3)
CR1<-brewer.pal(n = 3, name = 'YlOrRd')
CR2<-brewer.pal(n = 3, name = 'Blues')
Mcols<-c(CR1[parm.k.1],CR2[parm.k.2])

Mpch <- c(rep(21,4),rep(22,4))
Mpch2 <- c(NA,NA,NA,21,NA,NA,22,NA) # overlay symbols to differentiate equal k models

minAICcs <- apply(AICcs, 1, min)
delV.AICc <- AICcs - minAICcs
rnks.AICc <- t(apply(delV.AICc, 1, rank, ties.method='first'))
colnames(rnks.AICc) <- colnames(AICcs)

# Define delta AICc cut-off for "indistinguishably well performing" models
cut.AIC <- 2

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_ranks.pdf',
    height=6.5,width=2)
  par(mar=c(2,6.5,4,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
    plot(1:nrow(rnks.AICc), 1:nrow(rnks.AICc),
         type='n', yaxt='n',
         xlim=c(1,ncol(rnks.AICc)),
         ylim=c(0,nrow(rnks.AICc)+1),
         xlab='',
         ylab='',
         axes=F)
    title(xlab='Model rank by AICc',line=1)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
    axis(2, at=1:nrow(rnks.AICc), labels=labels, cex.axis=0.5, las=2)
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0))

    # Which models have delta-AICc within X=2 of best-performing model?
    xats <-table(which(delV.AICc < cut.AIC, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))+0.5
    segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
  
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.5
    polygon(pxats,pyats,col='grey90',border=NA)
    
    for(m in 1:ncol(rnks.AICc)){
      points(rnks.AICc[,m], 1:nrow(rnks.AICc), 
             type='p',  col='black', 
             bg=Mcols[m], pch=Mpch[m],
             cex=1, lwd=0.2)
      # Overlay secondary symbos
      points(rnks.AICc[,m], 1:nrow(rnks.AICc), 
             type='p',  col='black', 
             bg='white',pch=Mpch2[m],
             cex=0.3, lwd=0.2)
    }  
    box(lwd=1)
    
  par(xpd=TRUE)
  vadj <- 7.5
  xadj <- 1.5
  mods<-c(1,5)
    legend(xadj,nrow(rnks.AICc)+vadj,legend=colnames(rnks.AICc)[mods],
           pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
            pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8, 
           title=expression(italic(k)==1), bty='n')
    # overlay secondary symbols
    legend(xadj,nrow(rnks.AICc)+vadj,legend=colnames(rnks.AICc)[mods],
           pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
           pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
           title=expression(italic(k)==1), bty='n')
    
    mods<-c(2,6,7)
    legend(xadj+2,nrow(rnks.AICc)+vadj,legend=colnames(rnks.AICc)[mods],
           pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
           pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
           title=expression(italic(k)==2), bty='n')
    # overlay secondary symbols
    legend(xadj+2,nrow(rnks.AICc)+vadj,legend=colnames(rnks.AICc)[mods],
           pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
           pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
           title=expression(italic(k)==2), bty='n')
    
    mods<-c(3,4,8)
    legend(xadj+4,nrow(rnks.AICc)+vadj,legend=colnames(rnks.AICc)[mods],
           pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
           pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
           title=expression(italic(k)==3), bty='n')
    # overlay secondary symbols
    legend(xadj+4,nrow(rnks.AICc)+vadj,legend=colnames(rnks.AICc)[mods],
           pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
           pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
           title=expression(italic(k)==3), bty='n')
    
    rect(1.5,nrow(rnks.AICc)+vadj+1.75,7.5,nrow(rnks.AICc)+vadj-6)
    text(xadj+3,nrow(rnks.AICc)+vadj+0.75,'Models', adj=0.5, cex=0.8)
dev.off()

# ~~~~~~~~~
# Count times each model is in each rank
Cnt<-apply(rnks.AICc,2,function(x){table(factor(x,levels=1:ncol(rnks.AICc)))})
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

latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings.tex',label='table:AICc_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

# latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings.tex',label='table:AICc_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency) for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

setwd(wd)
# ~~~~~~~~~

# Of the times that AA is ranked first, how often are BD and CM within X AICc units?
cnt<-apply(delV.AICc[rnks.AICc[,ncol(rnks.AICc)]==1,3:4] < cut.AIC, 2, sum)
cnt
Cnt[1,'AA']
cnt/Cnt[1,'AA']
# Concl: BD is within 2 dAICc units > 1/3 the time, and CM is > 1/4 of the time

# Of times that AA is ranked first, how often are *either* BD and CM within X AICc units?
cnt<-sum(apply(delV.AICc[rnks.AICc[,ncol(rnks.AICc)]==1,3:4] < cut.AIC, 1, sum)>0)
cnt
cnt/Cnt[1,'AA']
# Concl: BD or CM are within 2 dAICc units ~60% of the time

# How many times is a single model the only best model?
cnt<-sum(apply(delV.AICc < cut.AIC, 1, sum)==1)
cnt/nrow(delV.AICc)

# ~~~~~~~~~
# What about for datasets that have a sample size of at least X?
SScuts <- seq(5,300,by=1)
fFirst.AICc<-fSecnd.AICc<-dim(0)
for(SScut in SScuts){
  Cnt<-apply(rnks.AICc[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnks.AICc)))})
  pCnt <- prop.table(Cnt,2)
  fFirst.AICc <- rbind(fFirst.AICc,pCnt[1,])
  fSecnd.AICc <- rbind(fSecnd.AICc,pCnt[2,])
}

# Darken colors of linear models
Mcols[c(1,5)]<-'grey40'

ltys <- c(1,1,1,6,2,2,3,1)
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_rankBySS.pdf',
    height=4,width=3)
par(mar=c(2,6.5,4,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.8)
  layout(matrix(c(1,2,3),ncol=1), heights=c(1,3,3)) 
  par(mar=c(0,3,0,0))
  plot(1,1,type="n", axes=F, xlab="", ylab="")
  
  xadj<-0.7
  yadj<-1.2
  lcex<-0.7
  
  mods<-c(1,5)
  legend(xadj,yadj,legend=colnames(fFirst.AICc)[mods],
         lty=ltys[mods], col=Mcols[mods],lwd=1.5,
         cex=lcex, ncol=1, pt.lwd=0.8, 
         title=expression(italic(k)==1), bty='n')
  
  mods<-c(2,6,7)
  legend(xadj+0.2,yadj,legend=colnames(fFirst.AICc)[mods],
         lty=ltys[mods], col=Mcols[mods],lwd=1.5,
         pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
         title=expression(italic(k)==2), bty='n')
  
  mods<-c(3,4,8)
  legend(xadj+0.4,yadj,legend=colnames(fFirst.AICc)[mods],
         lty=ltys[mods], col=Mcols[mods],lwd=1.5,
         pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
         title=expression(italic(k)==3), bty='n')
  
  text(xadj+0.3, yadj+0.1, 'Models', adj=0.5, cex=1)
  
  
  par(mar=c(3,3,0.5,1))
  matplot(SScuts,fFirst.AICc,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
          ylim=c(0,0.6),lwd=1.5,
          ylab='Fraction in first rank',
          xlab='Sample size greater than...')
  legend('topleft',legend='A',inset=-0.05,bty='n',cex=1.2)
  matplot(SScuts,fSecnd.AICc,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
          ylim=c(0,0.8),lwd=1.5,
          ylab='Fraction in second rank',
          xlab='Sample size greater than...')
  legend('topleft',legend='B',inset=-0.05,bty='n',cex=1.2)
dev.off()



SScut <- 50
Cnt<-apply(rnks.AICc[sample.sizes>=SScut,],2,function(x){table(factor(x,levels=1:ncol(rnks.AICc)))})
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

latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings_top50.tex',label='table:AICc_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

# latex(tab_Cnt,file='OnePredOnePrey_AICc_rankings_top50.tex',label='table:AICc_rankings_top50', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets (and percentage frequency)  with a sample size greater than 50 for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')

setwd(wd)
# ~~~~~~~~~
