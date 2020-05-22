# This script combines what used to be two seperate scripts
# ("plot_AICc.R" and "plot_RMSD.R") into a single file in order to 
# produce figures that contain both.

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
# Grab summary of AICc and RMSE estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

# For AICc
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

# Repeat for RMSD
RMSD.H1 <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Holling.I'][[1]][stat]}))
RMSD.H2 <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Holling.II'][[1]][stat]}))
RMSD.BD <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Beddington.DeAngelis'][[1]][stat]}))
RMSD.CM <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Crowley.Martin'][[1]][stat]}))
RMSD.R <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Ratio'][[1]][stat]}))
RMSD.AG <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Arditi.Ginzburg'][[1]][stat]}))
RMSD.HV <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Hassell.Varley'][[1]][stat]}))
RMSD.AA <- unlist(lapply(ffr.fits, function(x){ x$RMSD['Arditi.Akcakaya'][[1]][stat]}))

RMSDs <- data.frame(RMSD.H1, RMSD.H2, RMSD.BD, RMSD.CM, RMSD.R, RMSD.AG, RMSD.HV, RMSD.AA)
colnames(RMSDs) <- sub('RMSD.', '', colnames(RMSDs))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-AICc and ranks
minAICcs <- apply(AICcs, 1, min)
delV.AICc <- AICcs - minAICcs
rnks.AICc <- t(apply(delV.AICc, 1, rank, ties.method='first'))
colnames(rnks.AICc) <- colnames(AICcs)

# Define delta AICc cut-off for "indistinguishably well performing" models
cut.AIC <- 2

# "Delta-RMSD" and ranks
minRMSDs <- apply(RMSDs, 1, min)
delV.RMSD <- RMSDs - minRMSDs
rnks.RMSD <- t(apply(delV.RMSD, 1, rank, ties.method='first'))
colnames(rnks.RMSD) <- colnames(RMSDs)

# Define delta RMSD cut-off for "well performing" models
# Which models have an RMSD that is less than 1% of the data average (repeat for raw and bootstrapped datasets, which will each throw errors where the other does not apply, and merge)
cut.RMSD <- unlist(lapply(ffr.fits, 
                          function(x){0.01*mean(x$study.info$data.Nconsumed)}))
cut.RMSD2 <- unlist(lapply(ffr.fits, 
                           function(x){0.01*mean(x$study.info$data.Nconsumed.mean)}))
cut.RMSD[which(is.na(cut.RMSD))] <- cut.RMSD2[which(!is.na(cut.RMSD2))]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~
# Figure of Rank orders
#~~~~~~~~~~~~~~~~~~~~~~
# Define visual attributes for each model
parm.k.1<-c(1,2,3,3) # for Holling models
parm.k.2<-c(1,2,2,3) # for ratio models
CR1<-brewer.pal(n = 3, name = 'YlOrRd')
CR2<-brewer.pal(n = 3, name = 'Blues')
Mcols<-c(CR1[parm.k.1],CR2[parm.k.2])
Mpch <- c(rep(21,4),rep(22,4))
Mpch2 <- c(NA,NA,NA,21,NA,NA,22,NA) # overlay symbols to differentiate equal k models

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_and_RMSD_ranks.pdf',
    height=6.5,width=3)
nf<-layout(matrix(c(1,1,2,3),ncol=2,byrow=T), heights=c(0.7,6.3,6.3), widths=c(2,2,2))
# layout.show(nf)

# Legend ~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(0,0,0,0))
plot(0,0,ann=F,type='n',axes=F)

  vadj <- 0.7
  xadj1 <- -0.3
  xadj2 <- -0.1
  xadj3 <- 0.1
  
  mods<-c(1,5)
  legend(xadj1,vadj,legend=colnames(rnks.AICc)[mods],
         pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
         pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8, 
         title=expression(italic(k)==1), bty='n')
  # overlay secondary symbols
  legend(xadj1,vadj,legend=colnames(rnks.AICc)[mods],
         pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
         pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
         title=expression(italic(k)==1), bty='n')
  
  mods<-c(2,6,7)
  legend(xadj2,vadj,legend=colnames(rnks.AICc)[mods],
         pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
         pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
         title=expression(italic(k)==2), bty='n')
  # overlay secondary symbols
  legend(xadj2,vadj,legend=colnames(rnks.AICc)[mods],
         pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
         pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
         title=expression(italic(k)==2), bty='n')
  
  mods<-c(3,4,8)
  legend(xadj3,vadj,legend=colnames(rnks.AICc)[mods],
         pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
         pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
         title=expression(italic(k)==3), bty='n')
  # overlay secondary symbols
  legend(xadj3,vadj,legend=colnames(rnks.AICc)[mods],
         pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
         pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
         title=expression(italic(k)==3), bty='n')
  
  rect(xadj2-0.23,vadj+0.3,xadj2+0.43,vadj-1.55)
  text(xadj2,vadj+0.1,'Models', adj=0, cex=0.8)
  
  text(-0.7,-0.9,'(A) AICc',cex=0.75)
  text(0.7,-0.9,'(B) RMSD',cex=0.75)


# AICc results ~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(2,0.1,0,3.2))
par(mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
plot(1:nrow(rnks.AICc), 1:nrow(rnks.AICc),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnks.AICc)),
       ylim=c(0,nrow(rnks.AICc)+1),
       xlab='',
       ylab='',
       axes=F)
  title(xlab='Model rank',line=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
  axis(4, at=1:nrow(rnks.AICc), labels=labels, cex.axis=0.5, las=2, 
       hadj=0.5, mgp=c(0,3.2,0))
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


# RMSD results ~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(2,3.2,0,0.1))
par(mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
plot(1:nrow(rnks.RMSD), 1:nrow(rnks.RMSD),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnks.RMSD)),
       ylim=c(0,nrow(rnks.RMSD)+1),
       xlab='',
       ylab='',
       axes=F)
  title(xlab='Model rank',line=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white")
  axis(2, at=1:nrow(rnks.RMSD), labels=NA, cex.axis=0.5, las=2)
  axis(1, cex.axis=0.7, mgp=c(1.25,0,0),las=1)
  
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

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overall summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Count times each model is in each rank
Cnt_AICc<-apply(rnks.AICc,2,function(x){table(factor(x,levels=1:ncol(rnks.AICc)))})
Cnt_AICc
pCnt_AICc <- round(prop.table(Cnt_AICc,2)*100,1)
pCnt_AICc

Cnt_RMSD<-apply(rnks.RMSD,2,function(x){table(factor(x,levels=1:ncol(rnks.RMSD)))})
Cnt_RMSD
pCnt_RMSD <- round(prop.table(Cnt_RMSD,2)*100,1)
pCnt_RMSD

# Either pass counts or combine counts and proportions before exporting table
tab_Cnt_RMSD <- Cnt_RMSD
# bCnt_RMSD <- paste0(Cnt_RMSD,' (',pCnt_RMSD,')')
# bCnt_RMSD[Cnt_RMSD==0] <- 0
# tab_Cnt_RMSD <- matrix(bCnt_RMSD, nrow=nrow(Cnt_RMSD), byrow=T, dimnames=dimnames(Cnt_RMSD))

tab_Cnt_AICc <- Cnt_AICc
# bCnt_AICc <- paste0(Cnt_AICc,' (',pCnt_AICc,')')
# bCnt_AICc[Cnt_AICc==0] <- 0
# tab_Cnt_AICc <- matrix(bCnt_AICc, nrow=nrow(Cnt_AICc), byrow=T, dimnames=dimnames(Cnt_AICc))

tab_Cnt <- rbind(tab_Cnt_AICc, rep('',ncol(tab_Cnt_AICc)), tab_Cnt_RMSD)

# ~~~~~~~~~~~~~~~
# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')

latex(tab_Cnt,file='OnePredOnePrey_AICc_and_RMSD_rankings.tex',label='table:AICc_and_RMSD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$ or by RMSD.')

# latex(tab_Cnt_AICc,file='OnePredOnePrey_AICc_rankings.tex',label='table:AICc_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')
# 
# latex(tab_Cnt_RMSD,file='OnePredOnePrey_RMSD_rankings.tex',label='table:RMSD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by the root mean square deviation (RMSD).')

setwd(wd)

# ~~~~~~~~~~~~~~~~~~~~~~~
# More summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~
# Of the times that AA is ranked first, how often are BD and CM within 2 AICc units?
cnt_AICc<-apply(delV.AICc[rnks.AICc[,ncol(rnks.AICc)]==1,3:4] < cut.AIC, 2, sum)
cnt_AICc
Cnt_AICc[1,'AA']
cnt_AICc/Cnt_AICc[1,'AA']

cnt_RMSD<-apply(delV.RMSD[rnks.RMSD[,ncol(rnks.RMSD)]==1,3:4] < cut.RMSD, 2, sum)
cnt_RMSD
Cnt_RMSD[1,'AA']
cnt_RMSD/Cnt_RMSD[1,'AA']

# Of times that AA is ranked first, how often are *either* BD and CM within 2 AICc units?
cnt_AICc<-sum(apply(delV.AICc[rnks.AICc[,ncol(rnks.AICc)]==1,3:4] < cut.AIC, 1, sum)>0)
cnt_AICc
cnt_AICc/Cnt_AICc[1,'AA']

cnt_RMSD<-sum(apply(delV.RMSD[rnks.RMSD[,ncol(rnks.RMSD)]==1,3:4] < cut.RMSD, 1, sum)>0)
cnt_RMSD
cnt_RMSD/Cnt_RMSD[1,'AA']

# How many times is a single model the only best model?
cnt_AICc<-sum(apply(delV.AICc < cut.AIC, 1, sum)==1)
cnt_AICc/nrow(delV.AICc)

# ~~~~~~~~~
# What about for datasets that have a sample size of at least X?
SScuts <- seq(5,300,by=1)
fFirst.AICc<-fSecnd.AICc<-fFirst.RMSD<-fSecnd.RMSD<-dim(0)
for(SScut in SScuts){
  Cnt_AICc<-apply(rnks.AICc[sample.sizes>=SScut,],2, 
                  function(x){
                    table(factor(x,levels=1:ncol(rnks.AICc)))})
  pCnt_AICc <- prop.table(Cnt_AICc,2)
  Cnt_RMSD<-apply(rnks.RMSD[sample.sizes>=SScut,],2, 
                  function(x){
                    table(factor(x,levels=1:ncol(rnks.RMSD)))})
  pCnt_RMSD <- prop.table(Cnt_RMSD,2)
  
  fFirst.AICc <- rbind(fFirst.AICc,pCnt_AICc[1,])
  fSecnd.AICc <- rbind(fSecnd.AICc,pCnt_AICc[2,])
  fFirst.RMSD <- rbind(fFirst.RMSD,pCnt_RMSD[1,])
  fSecnd.RMSD <- rbind(fSecnd.RMSD,pCnt_RMSD[2,])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure of top two rankings as a function of sample size
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Darken colors of linear models relative to above
Mcols[c(1,5)]<-'grey40'
ltys <- c(1,1,1,6,2,2,3,1)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_and_RMSD_toprankBySS.pdf',
    height=4,width=5)
par(mar=c(2,6.5,4,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.8)
nf<-layout(matrix(c(1,1,2,4,3,5),ncol=2,byrow=T), heights=c(1.2,3,3)) 
# layout.show(nf)

# legend ~~~~~~~~~
par(mar=c(0,0,0,0))
plot(0,0,type="n", axes=F, xlab="", ylab="")

xadj1 <- -0.2
xadj2 <- 0
xadj3 <- 0.2
yadj <- 0.5
lcex <- 0.7

mods<-c(1,5)
legend(xadj1,yadj,legend=colnames(fFirst.AICc)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       cex=lcex, ncol=1, pt.lwd=0.8, 
       title=expression(italic(k)==1), bty='n')

mods<-c(2,6,7)
legend(xadj2,yadj,legend=colnames(fFirst.AICc)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==2), bty='n')

mods<-c(3,4,8)
legend(xadj3,yadj,legend=colnames(fFirst.AICc)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==3), bty='n')

rect(xadj2-0.2,vadj+0.25,xadj2+0.4,vadj-1.7)
text(xadj2+0.02, yadj+0.18, 'Models', adj=0,cex=1)

par(mar=c(2.5,3,0.5,0))
matplot(SScuts,fFirst.AICc,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
        ylim=c(0,0.7),lwd=1.5,
        ylab='Fraction in first rank',
        xlab='')
  mtext('A',side=3,line=-1.25,at=5,cex=0.9)
  mtext('AICc',side=3,line=0.3,at=150,cex=1)
  
par(mar=c(3,3,0,0))
matplot(SScuts,fSecnd.AICc,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
        ylim=c(0,0.7),lwd=1.5,
        ylab='Fraction in second rank',
        xlab='Sample size greater than...')
  mtext('B',side=3,line=-1.25,at=5,cex=0.9)
  
par(mar=c(2.5,2,0.5,1))
matplot(SScuts,fFirst.RMSD,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
        ylim=c(0,0.7),lwd=1.5,
        ylab='',
        xlab='')
  mtext('C',side=3,line=-1.25,at=5,cex=0.9)
  mtext('RMSD',side=3,line=0.3,at=150,cex=1)
  
par(mar=c(3,2,0,1))
matplot(SScuts,fSecnd.RMSD,las=1,type='l', col=Mcols,lty=ltys,cex=0.5,
        ylim=c(0,0.7),lwd=1.5,
        ylab='',
        xlab='Sample size greater than...')
  mtext('D',side=3,line=-1.25,at=5,cex=0.9)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~