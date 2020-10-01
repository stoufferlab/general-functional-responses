source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

library(RColorBrewer)
library(stats) # for rug()
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
labels <- paste0(labels, ' (',sample.sizes,')')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grab summary of AIC and BIC estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

# For AIC
AIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.I'][[1]][stat]}))
AIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.II'][[1]][stat]}))
AIC.BD <- unlist(lapply(ffr.fits, function(x){ x$AIC['Beddington.DeAngelis'][[1]][stat]}))
AIC.CM <- unlist(lapply(ffr.fits, function(x){ x$AIC['Crowley.Martin'][[1]][stat]}))
AIC.R <- unlist(lapply(ffr.fits, function(x){ x$AIC['Ratio'][[1]][stat]}))
AIC.HV <- unlist(lapply(ffr.fits, function(x){ x$AIC['Hassell.Varley'][[1]][stat]}))
AIC.AG <- unlist(lapply(ffr.fits, function(x){ x$AIC['Arditi.Ginzburg'][[1]][stat]}))
AIC.AA <- unlist(lapply(ffr.fits, function(x){ x$AIC['Arditi.Akcakaya'][[1]][stat]}))

AICs <- data.frame(AIC.H1, AIC.H2, AIC.BD, AIC.CM, AIC.R, AIC.HV, AIC.AG, AIC.AA)
colnames(AICs) <- sub('AIC.', '', colnames(AICs))
# Hack to rename R -> LR (ratio -> linear ratio)
colnames(AICs) <- sub('R','LR',colnames(AICs))

# For BIC
BIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$BIC['Holling.I'][[1]][stat]}))
BIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$BIC['Holling.II'][[1]][stat]}))
BIC.BD <- unlist(lapply(ffr.fits, function(x){ x$BIC['Beddington.DeAngelis'][[1]][stat]}))
BIC.CM <- unlist(lapply(ffr.fits, function(x){ x$BIC['Crowley.Martin'][[1]][stat]}))
BIC.R <- unlist(lapply(ffr.fits, function(x){ x$BIC['Ratio'][[1]][stat]}))
BIC.HV <- unlist(lapply(ffr.fits, function(x){ x$BIC['Hassell.Varley'][[1]][stat]}))
BIC.AG <- unlist(lapply(ffr.fits, function(x){ x$BIC['Arditi.Ginzburg'][[1]][stat]}))
BIC.AA <- unlist(lapply(ffr.fits, function(x){ x$BIC['Arditi.Akcakaya'][[1]][stat]}))

BICs <- data.frame(BIC.H1, BIC.H2, BIC.BD, BIC.CM, BIC.R, BIC.HV, BIC.AG, BIC.AA)
colnames(BICs) <- sub('BIC.', '', colnames(BICs))
# Hack to rename R -> LR (ratio -> linear ratio)
colnames(BICs) <- sub('R','LR',colnames(BICs))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-AIC and ranks
minAICs <- apply(AICs, 1, min)
delta.AIC <- AICs - minAICs
rnks.AIC <- t(apply(AICs, 1, rank, ties.method='first'))
colnames(rnks.AIC) <- colnames(AICs)

# Define delta AIC cut-off for "indistinguishably well performing" models
cutoff.AIC <- 2

# Delta-BIC and ranks
minBICs <- apply(BICs, 1, min)
delta.BIC <- BICs - minBICs
rnks.BIC <- t(apply(BICs, 1, rank, ties.method='first'))
colnames(rnks.BIC) <- colnames(BICs)

# Define delta BIC cut-off for "indistinguishably well performing" models
cutoff.BIC <- 2

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

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AIC_and_BIC_ranks.pdf',
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
  legend(xadj1,vadj,legend=colnames(rnks.AIC)[mods],
         pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
         pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8, 
         title=expression(italic(k)==1), bty='n')
  # overlay secondary symbols
  legend(xadj1,vadj,legend=colnames(rnks.AIC)[mods],
         pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
         pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
         title=expression(italic(k)==1), bty='n')
  
  mods<-c(2,6,7)
  legend(xadj2,vadj,legend=colnames(rnks.AIC)[mods],
         pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
         pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
         title=expression(italic(k)==2), bty='n')
  # overlay secondary symbols
  legend(xadj2,vadj,legend=colnames(rnks.AIC)[mods],
         pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
         pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
         title=expression(italic(k)==2), bty='n')
  
  mods<-c(3,4,8)
  legend(xadj3,vadj,legend=colnames(rnks.AIC)[mods],
         pch=Mpch[mods], pt.bg=Mcols[mods], col='black', bg='white',
         pt.cex=1.1,cex=0.6, ncol=1, pt.lwd=0.8,
         title=expression(italic(k)==3), bty='n')
  # overlay secondary symbols
  legend(xadj3,vadj,legend=colnames(rnks.AIC)[mods],
         pch=Mpch2[mods], pt.bg='white', col='black', bg=NA,
         pt.cex=0.3,cex=0.6, ncol=1, pt.lwd=0.3,
         title=expression(italic(k)==3), bty='n')
  
  rect(xadj2-0.23,vadj+0.3,xadj2+0.43,vadj-1.55)
  text(xadj2,vadj+0.1,'Models', adj=0, cex=0.8)
  
  text(-0.7,-0.9,'(A) AIC',cex=0.75)
  text(0.7,-0.9,'(B) BIC',cex=0.75)


# AIC results ~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(2,0.1,0,3.2))
par(mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
plot(1:nrow(rnks.AIC), 1:nrow(rnks.AIC),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnks.AIC)),
       ylim=c(0,nrow(rnks.AIC)+1),
       xlab='',
       ylab='',
       axes=F)
  title(xlab='Model rank',line=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
  axis(4, at=1:nrow(rnks.AIC), labels=labels, cex.axis=0.5, las=2, 
       hadj=0.5, mgp=c(0,3.2,0))
  axis(1, cex.axis=0.7, mgp=c(1.25,0,0))
  
  # Which models have delta-AIC within X=2 of best-performing model?
  xats <-table(factor(which(delta.AIC < cutoff.AIC, arr.ind=T)[,1],
                      levels=1:nrow(delta.AIC)))+0.5
  yats <- 0:(length(xats))+0.5
  segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
  segments(c(0.5, xats[-length(xats)]), yats, 
           c(0.5, xats[-1]), yats, col='black')
  
  pxats<-c(0,rep(xats,each=2),0)
  pyats<-rep(0:(length(xats)),each=2)+0.5
  polygon(pxats,pyats,col='grey90',border=NA)
  
  for(m in 1:ncol(rnks.AIC)){
    points(rnks.AIC[,m], 1:nrow(rnks.AIC), 
           type='p',  col='black', 
           bg=Mcols[m], pch=Mpch[m],
           cex=1, lwd=0.2)
    # Overlay secondary symbos
    points(rnks.AIC[,m], 1:nrow(rnks.AIC), 
           type='p',  col='black', 
           bg='white',pch=Mpch2[m],
           cex=0.3, lwd=0.2)
  }  
  box(lwd=1)


# BIC results ~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(2,3.2,0,0.1))
par(mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
plot(1:nrow(rnks.BIC), 1:nrow(rnks.BIC),
     type='n', yaxt='n',
     xlim=c(1,ncol(rnks.BIC)),
     ylim=c(0,nrow(rnks.BIC)+1),
     xlab='',
     ylab='',
     axes=F)
  title(xlab='Model rank',line=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
  axis(4, at=1:nrow(rnks.BIC), labels=labels, cex.axis=0.5, las=2, 
       hadj=0.5, mgp=c(0,3.2,0))
  axis(1, cex.axis=0.7, mgp=c(1.25,0,0))

  # Which models have delta-BIC within X=2 of best-performing model?
  xats <-table(factor(which(delta.BIC < cutoff.BIC, arr.ind=T)[,1],
                      levels=1:nrow(delta.BIC)))+0.5
  yats <- 0:(length(xats))+0.5
  segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
  segments(c(0.5, xats[-length(xats)]), yats, 
           c(0.5, xats[-1]), yats, col='black')
  
  pxats<-c(0,rep(xats,each=2),0)
  pyats<-rep(0:(length(xats)),each=2)+0.5
  polygon(pxats,pyats,col='grey90',border=NA)
  
  for(m in 1:ncol(rnks.BIC)){
    points(rnks.BIC[,m], 1:nrow(rnks.BIC), 
           type='p',  col='black', 
           bg=Mcols[m], pch=Mpch[m],
           cex=1, lwd=0.2)
    # Overlay secondary symbos
    points(rnks.BIC[,m], 1:nrow(rnks.BIC), 
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
mod.order <- colnames(rnks.AIC)
mods.n <- length(mod.order)
Cnt_AIC <- apply(rnks.AIC,2,function(x){table(factor(x,levels=1:mods.n))})
Cnt_AIC <- Cnt_AIC[,mod.order]
Cnt_AIC
sum(Cnt_AIC[1,c('BD','CM','AA')])
pCnt_AIC <- round(prop.table(Cnt_AIC,2)*100,1)
pCnt_AIC

mod.order <- colnames(rnks.BIC)
mods.n <- length(mod.order)
Cnt_BIC <- apply(rnks.BIC,2,function(x){table(factor(x,levels=1:mods.n))})
Cnt_BIC <- Cnt_BIC[,mod.order]
Cnt_BIC
sum(Cnt_BIC[1,c('BD','CM','AA')])
pCnt_BIC <- round(prop.table(Cnt_BIC,2)*100,1)
pCnt_BIC

# Either pass counts or combine counts and proportions before exporting table
tab_Cnt_AIC <- Cnt_AIC
# bCnt_AIC <- paste0(Cnt_AIC,' (',pCnt_AIC,')')
# bCnt_AIC[Cnt_AIC==0] <- 0
# tab_Cnt_AIC <- matrix(bCnt_AIC, nrow=nrow(Cnt_AIC), byrow=T, dimnames=dimnames(Cnt_AIC))

tab_Cnt_BIC <- Cnt_BIC
# bCnt_BIC <- paste0(Cnt_BIC,' (',pCnt_BIC,')')
# bCnt_BIC[Cnt_BIC==0] <- 0
# tab_Cnt_BIC <- matrix(bCnt_BIC, nrow=nrow(Cnt_BIC), byrow=T, dimnames=dimnames(Cnt_BIC))


tab_Cnt <- rbind(c(rep('',3),'$AIC$',rep('',4)),
                 tab_Cnt_AIC, 
                 rep('',ncol(tab_Cnt_AIC)), 
                 c(rep('',3),'$BIC$',rep('',4)),
                 tab_Cnt_BIC)

# ~~~~~~~~~~~~~~~
# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')

latex(tab_Cnt,file='OnePredOnePrey_AIC_and_BIC_rankings.tex',label='table:AIC_and_BIC_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC$ or by $BIC$.')

setwd(wd)

# ~~~~~~~~~~~~~~~~~~~~~~~
# More summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~
# How many times is a single model the only best model by AIC?
cnt_AIC_single<-sum(apply(delta.AIC <= cutoff.AIC, 1, sum)==1)
cnt_AIC_single
cnt_AIC_single/nrow(delta.AIC)*100

#~~~~~
# How often are *each of* the other models within 2 AIC units of the top model?
r1Mod.AIC <- which.max(Cnt_AIC[1,])
cnt_AIC<-apply(delta.AIC[rnks.AIC[,r1Mod.AIC]==1,-r1Mod.AIC] < cutoff.AIC, 2, sum)
cnt_AIC

# How often is any other model within 2 AIC units of the top model?
cnt_AIC<-sum(apply(delta.AIC[rnks.AIC[,r1Mod.AIC]==1,-r1Mod.AIC] < cutoff.AIC, 1, sum)>0)
cnt_AIC
cnt_AIC/Cnt_AIC[1,r1Mod.AIC]*100

#~~~~~
# How often are other models within 2 AIC units of the 2nd-best model?
r2Mod.AIC <- which(rank(-Cnt_AIC[1,])==2)
cnt_AIC<-apply(delta.AIC[rnks.AIC[,r2Mod.AIC]==1,-r2Mod.AIC] < cutoff.AIC, 2, sum)
cnt_AIC

# How often is any other model within 2 AIC units of the 2nd-best model?
cnt_AIC<-sum(apply(delta.AIC[rnks.AIC[,r2Mod.AIC]==1,-r2Mod.AIC] < cutoff.AIC, 1, sum)>0)
cnt_AIC
cnt_AIC/Cnt_AIC[1,r2Mod.AIC]

#~~~~~
# How often are other models within 2 AIC units of the 3rd-best model?
r3Mod.AIC <- which(rank(-Cnt_AIC[1,])==3)
cnt_AIC<-apply(delta.AIC[rnks.AIC[,r3Mod.AIC]==1,-r3Mod.AIC] < cutoff.AIC, 2, sum)
cnt_AIC

# How often is any other model within 2 AIC units of the 3rd-best model?
cnt_AIC<-sum(apply(delta.AIC[rnks.AIC[,r3Mod.AIC]==1,-r3Mod.AIC] < cutoff.AIC, 1, sum)>0)
cnt_AIC
cnt_AIC/Cnt_AIC[1,r3Mod.AIC]

#~~~~~~~~~~~~~~~~~~~~~~

# How many times is a single model the only best model by BIC?
cnt_BIC_single<-sum(apply(delta.BIC <= cutoff.BIC, 1, sum)==1)
cnt_BIC_single
cnt_BIC_single/nrow(delta.BIC)*100

#~~~~~
# How often are *each of* the other models within 2 BIC units of the top model?
r1Mod.BIC <- which.max(Cnt_BIC[1,])
cnt_BIC<-apply(delta.BIC[rnks.BIC[,r1Mod.BIC]==1,-r1Mod.BIC] < cutoff.BIC, 2, sum)
cnt_BIC

# How often is any other model within 2 BIC units of the top model?
cnt_BIC<-sum(apply(delta.BIC[rnks.BIC[,r1Mod.BIC]==1,-r1Mod.BIC] < cutoff.BIC, 1, sum)>0)
cnt_BIC
cnt_BIC/Cnt_BIC[1,r1Mod.BIC]*100

#~~~~~
# How often are other models within 2 BIC units of the 2nd-best model?
r2Mod.BIC <- which(rank(-Cnt_BIC[1,])==2)
cnt_BIC<-apply(delta.BIC[rnks.BIC[,r2Mod.BIC]==1,-r2Mod.BIC] < cutoff.BIC, 2, sum)
cnt_BIC

# How often is any other model within 2 BIC units of the 2nd-best model?
cnt_BIC<-sum(apply(delta.BIC[rnks.BIC[,r2Mod.BIC]==1,-r2Mod.BIC] < cutoff.BIC, 1, sum)>0)
cnt_BIC
cnt_BIC/Cnt_BIC[1,r2Mod.BIC]

#~~~~~
# How often are other models within 2 BIC units of the 3rd-best model?
r3Mod.BIC <- which(rank(-Cnt_BIC[1,])==3)
cnt_BIC<-apply(delta.BIC[rnks.BIC[,r3Mod.BIC]==1,-r3Mod.BIC] < cutoff.BIC, 2, sum)
cnt_BIC

# How often is any other model within 2 BIC units of the 3rd-best model?
cnt_BIC<-sum(apply(delta.BIC[rnks.BIC[,r3Mod.BIC]==1,-r3Mod.BIC] < cutoff.BIC, 1, sum)>0)
cnt_BIC
cnt_BIC/Cnt_BIC[1,r3Mod.BIC]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What about for datasets that have a sample size of at least X?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# hist(sample.sizes,breaks=30)
# hist(log(sample.sizes),breaks=30)
range(sample.sizes)
mean(sample.sizes)
median(sample.sizes)

# ~~~~~~~~~~~~~~~~~~~

# when excluding smaller than
# mxSS <- 300
# mnSS <- min(sample.sizes)

# when excluding greater than
mxSS <- max(sample.sizes)
mnSS <- 20

SScuts <- seq(mnSS,mxSS,by=1)
fFirst.AIC<-fSecnd.AIC<-fFirst.BIC<-fSecnd.BIC<-dim(0)
pfFirst.AIC<-pfSecnd.AIC<-pfFirst.BIC<-pfSecnd.BIC<-dim(0)
for(SScut in SScuts){
  
  # Cnt_AIC<-apply(rnks.AIC[sample.sizes>=SScut,],2, # exclude smaller than
  #                 function(x){
  #                   table(factor(x,levels=1:mods.n))})
  Cnt_AIC<-apply(rnks.AIC[sample.sizes<=SScut,],2, # exclude greater than
                  function(x){
                    table(factor(x,levels=1:mods.n))})
  Cnt_AIC <- Cnt_AIC[,mod.order]
  pCnt_AIC <- prop.table(Cnt_AIC,2)
  
  # Cnt_BIC<-apply(rnks.BIC[sample.sizes>=SScut,],2, # exclude smaller than
  #                 function(x){
  #                   table(factor(x,levels=1:mods.n))})
  Cnt_BIC<-apply(rnks.BIC[sample.sizes<=SScut,],2, # exclude greater than
                 function(x){
                   table(factor(x,levels=1:mods.n))})
  Cnt_BIC <- Cnt_BIC[,mod.order]
  pCnt_BIC <- prop.table(Cnt_BIC,2)
  
  fFirst.AIC <- rbind(fFirst.AIC,Cnt_AIC[1,])
  fSecnd.AIC <- rbind(fSecnd.AIC,Cnt_AIC[2,])
  fFirst.BIC <- rbind(fFirst.BIC,Cnt_BIC[1,])
  fSecnd.BIC <- rbind(fSecnd.BIC,Cnt_BIC[2,])

  
  pfFirst.AIC <- rbind(pfFirst.AIC,pCnt_AIC[1,])
  pfSecnd.AIC <- rbind(pfSecnd.AIC,pCnt_AIC[2,])
  pfFirst.BIC <- rbind(pfFirst.BIC,pCnt_BIC[1,])
  pfSecnd.BIC <- rbind(pfSecnd.BIC,pCnt_BIC[2,])
  
}

# How many datasets are included at the mxSS=300 limit of the plots?
sum(sample.sizes>=mxSS)

# What are stats at the median sample size?
s=median(sample.sizes)
length(sample.sizes[sample.sizes<=s])

fFirst.AIC[which(SScuts==s),]
round(pfFirst.AIC[which(SScuts==s),],3)*100
fSecnd.AIC[which(SScuts==s),]
round(pfSecnd.AIC[which(SScuts==s),],3)*100

fFirst.BIC[which(SScuts==s),]
round(pfFirst.BIC[which(SScuts==s),],3)*100
fSecnd.BIC[which(SScuts==s),]
round(pfSecnd.BIC[which(SScuts==s),],3)*100

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure of top two rankings as a function of sample size
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Darken colors of linear models relative to above
Mcols[c(1,5)]<-'grey40'
ltys <- c(1,1,1,6,2,2,1,1)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AIC_and_BIC_toprankBySS.pdf',
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

ylims <- c(0,0.7)

mods<-c(1,5)
legend(xadj1,yadj,legend=colnames(pfFirst.AIC)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       cex=lcex, ncol=1, pt.lwd=0.8, 
       title=expression(italic(k)==1), bty='n')

mods<-c(2,6,7)
legend(xadj2,yadj,legend=colnames(pfFirst.AIC)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==2), bty='n')

mods<-c(3,4,8)
legend(xadj3,yadj,legend=colnames(pfFirst.AIC)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==3), bty='n')

rect(xadj2-0.2,yadj+0.4,xadj2+0.4,yadj-1.5)
text(xadj2+0.02, yadj+0.18, 'Models', adj=0,cex=1)

par(mar=c(2.5,3,0.5,0))
matplot(SScuts,pfFirst.AIC,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='Fraction in first rank',
        xlab='')
  rug(sample.sizes,side=3,quiet=TRUE)
  mtext('A',side=3,line=-1.5,at=mnSS,cex=0.9)
  mtext('AIC',side=3,line=0.3,adj=0.5,cex=1)

par(mar=c(3,3,0,0))
matplot(SScuts,pfSecnd.AIC,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='Fraction in second rank',
        # xlab='Sample size greater than...')
        xlab='Sample size less than...')
  rug(sample.sizes,side=3,quiet=TRUE)
  mtext('B',side=3,line=-1.5,at=mnSS,cex=0.9)

par(mar=c(2.5,2,0.5,1))
matplot(SScuts,pfFirst.BIC,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='',
        xlab='')
  rug(sample.sizes,side=3,quiet=TRUE)
  mtext('C',side=3,line=-1.5,at=mnSS,cex=0.9)
  mtext('BIC',side=3,line=0.3,adj=0.5,cex=1)

par(mar=c(3,2,0,1))
matplot(SScuts,pfSecnd.BIC,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='',
        # xlab='Sample size greater than...')
        xlab='Sample size less than...')
  rug(sample.sizes,side=3,quiet=TRUE)
  mtext('D',side=3,line=-1.5,at=mnSS,cex=0.9)

dev.off()

