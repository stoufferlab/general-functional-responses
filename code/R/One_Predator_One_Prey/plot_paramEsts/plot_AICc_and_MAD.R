source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

library(RColorBrewer)
library(stats) # for rug()
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
# Grab summary of AICc and MAD estimates across bootstrapped fits
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

# Repeat for MAD
MAD.H1 <- unlist(lapply(ffr.fits, function(x){ x$MAD['Holling.I'][[1]][stat]}))
MAD.H2 <- unlist(lapply(ffr.fits, function(x){ x$MAD['Holling.II'][[1]][stat]}))
MAD.BD <- unlist(lapply(ffr.fits, function(x){ x$MAD['Beddington.DeAngelis'][[1]][stat]}))
MAD.CM <- unlist(lapply(ffr.fits, function(x){ x$MAD['Crowley.Martin'][[1]][stat]}))
MAD.R <- unlist(lapply(ffr.fits, function(x){ x$MAD['Ratio'][[1]][stat]}))
MAD.AG <- unlist(lapply(ffr.fits, function(x){ x$MAD['Arditi.Ginzburg'][[1]][stat]}))
MAD.HV <- unlist(lapply(ffr.fits, function(x){ x$MAD['Hassell.Varley'][[1]][stat]}))
MAD.AA <- unlist(lapply(ffr.fits, function(x){ x$MAD['Arditi.Akcakaya'][[1]][stat]}))

MADs <- data.frame(MAD.H1, MAD.H2, MAD.BD, MAD.CM, MAD.R, MAD.AG, MAD.HV, MAD.AA)
colnames(MADs) <- sub('MAD.', '', colnames(MADs))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-AICc and ranks
minAICcs <- apply(AICcs, 1, min)
delta.AICc <- AICcs - minAICcs
rnks.AICc <- t(apply(AICcs, 1, rank, ties.method='first'))
colnames(rnks.AICc) <- colnames(AICcs)

# Define delta AICc cut-off for "indistinguishably well performing" models
cutoff.AIC <- 2

# "Delta-MAD" and ranks
minMADs <- apply(MADs, 1, min)
delta.MAD <- MADs - minMADs
rnks.MAD <- t(apply(MADs, 1, rank, ties.method='first'))
colnames(rnks.MAD) <- colnames(MADs)


# Define MAD cut-off for "well performing" models
# Calculate feeding count averaged across all treatment in order to later determine which models have an MAD that is less than x% of the data average (repeat for raw and bootstrapped datasets, which will each throw errors where the other does not apply, and then merge)
mean.Nconsumed <- unlist(lapply(ffr.fits, 
                          function(x){mean(x$study.info$data.Nconsumed)}))
mean.Nconsumed_boot <- unlist(lapply(ffr.fits, 
                           function(x){mean(x$study.info$data.Nconsumed.mean*x$study.info$data.n)}))
mean.Nconsumed[which(is.na(mean.Nconsumed))] <- mean.Nconsumed_boot[which(!is.na(mean.Nconsumed_boot))]

perc.cut <- 10
cutoff.MAD <- (perc.cut/100)*mean.Nconsumed

# Define delta-MAD cut-off ****relative on best-performing model****
perc.cut <- 1
cutoff.delta.MAD <- (perc.cut/100)*minMADs
  
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

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_and_MAD_ranks.pdf',
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
  text(0.7,-0.9,'(B) MAD',cex=0.75)


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
  xats <-table(factor(which(delta.AICc < cutoff.AIC, arr.ind=T)[,1],
                      levels=1:nrow(delta.AICc)))+0.5
  yats <- 0:(length(xats))+0.5
  segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
  segments(c(0.5, xats[-length(xats)]), yats, 
           c(0.5, xats[-1]), yats, col='black')
  
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


# MAD results ~~~~~~~~~~~~~~~~~~~~~~~~
par(mar=c(2,3.2,0,0.1))
par(mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
plot(1:nrow(rnks.MAD), 1:nrow(rnks.MAD),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnks.MAD)),
       ylim=c(0,nrow(rnks.MAD)+1),
       xlab='',
       ylab='',
       axes=F)
  title(xlab='Model rank',line=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white")
  axis(2, at=1:nrow(rnks.MAD), labels=NA, cex.axis=0.5, las=2)
  axis(1, cex.axis=0.7, mgp=c(1.25,0,0),las=1)
  
  # Which models have a MAD below the cut-off?
  # xats <-table(factor(which(MADs < cutoff.MAD, arr.ind=T)[,1],
  #              levels=1:nrow(MADs)))+0.5
  
  # Which models have a delta-MAD below the delta cut-off?
  xats <-table(factor(which(delta.MAD < cutoff.delta.MAD, arr.ind=T)[,1],
                      levels=1:nrow(MADs)))+0.5
  
  yats <- 0:(length(xats))+0.5
  segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
  segments(c(0.5, xats[-length(xats)]), yats, 
           c(0.5, xats[-1]), yats, col='black')
  
  pxats<-c(0,rep(xats,each=2),0)
  pyats<-rep(0:(length(xats)),each=2)+0.5
  polygon(pxats,pyats,col='grey90',border=NA)
  
  for(m in 1:ncol(rnks.MAD)){
    points(rnks.MAD[,m], 1:nrow(rnks.MAD), 
           type='p',  col='black', 
           bg=Mcols[m], pch=Mpch[m],
           cex=1, lwd=0.2)
    # Overlay secondary symbos
    points(rnks.MAD[,m], 1:nrow(rnks.MAD), 
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
# Total number of datasets
nrow(rnks.AICc)

# Count times each model is in each rank
mod.order <- colnames(rnks.AICc)
mods.n <- length(mod.order)
Cnt_AICc <- apply(rnks.AICc,2,function(x){table(factor(x,levels=1:mods.n))})
Cnt_AICc <- Cnt_AICc[,mod.order]
Cnt_AICc
sum(Cnt_AICc[1,c('BD','CM','AA')])
nrow(rnks.AICc)-sum(Cnt_AICc[1,c('BD','CM','AA')])
pCnt_AICc <- round(prop.table(Cnt_AICc,2)*100,1)
pCnt_AICc

Cnt_MAD <- apply(rnks.MAD,2,function(x){table(factor(x,levels=1:mods.n))})
Cnt_MAD <- Cnt_MAD[,mod.order]
Cnt_MAD
sum(Cnt_MAD[1,c('BD','CM','AA')])
pCnt_MAD <- round(prop.table(Cnt_MAD,2)*100,1)
pCnt_MAD

# Either pass counts or combine counts and proportions before exporting table
tab_Cnt_AICc <- Cnt_AICc
# bCnt_AICc <- paste0(Cnt_AICc,' (',pCnt_AICc,')')
# bCnt_AICc[Cnt_AICc==0] <- 0
# tab_Cnt_AICc <- matrix(bCnt_AICc, nrow=nrow(Cnt_AICc), byrow=T, dimnames=dimnames(Cnt_AICc))

tab_Cnt_MAD <- Cnt_MAD
# bCnt_MAD <- paste0(Cnt_MAD,' (',pCnt_MAD,')')
# bCnt_MAD[Cnt_MAD==0] <- 0
# tab_Cnt_MAD <- matrix(bCnt_MAD, nrow=nrow(Cnt_MAD), byrow=T, dimnames=dimnames(Cnt_MAD))

tab_Cnt <- rbind(c(rep('',3),'AICc',rep('',4)),
                 tab_Cnt_AICc, 
                 rep('',ncol(tab_Cnt_AICc)), 
                 c(rep('',3),'MAD',rep('',4)),
                 tab_Cnt_MAD)

# ~~~~~~~~~~~~~~~
# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')

latex(tab_Cnt,file='OnePredOnePrey_AICc_and_MAD_rankings.tex',label='table:AICc_and_MAD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$ or by $MAD$.')

# latex(tab_Cnt_AICc,file='OnePredOnePrey_AICc_rankings.tex',label='table:AICc_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by $AIC_c$.')
# 
# latex(tab_Cnt_MAD,file='OnePredOnePrey_MAD_rankings.tex',label='table:MAD_rankings', rowlabel='Rank', na.blank=TRUE, caption='The number of datasets for which each functional response model achieved a given rank relative to all other models as judged by the root mean square deviation (MAD).')

setwd(wd)

# ~~~~~~~~~~~~~~~~~~~~~~~
# More summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~
# How many times is a single model the only best model by AICc?
cnt_AICc_single<-sum(apply(delta.AICc <= cutoff.AIC, 1, sum)==1)
cnt_AICc_single
round(cnt_AICc_single/nrow(delta.AICc)*100,1)

#~~~~~
# How often are *each of* the other models within 2 AICc units of the top model?
Cnt_AICc[1,which(Cnt_AICc[1,]==max(Cnt_AICc[1,]))]
r1Mod.AICc <- 4 # for CM
# r1Mod.AICc <- 8 # for AA
cnt_AICc<-apply(delta.AICc[rnks.AICc[,r1Mod.AICc]==1,-r1Mod.AICc] < cutoff.AIC, 2, sum)
cnt_AICc

# How often is any other model within 2 AICc units of the top model?
r1Mod.AICc <- 4 # for CM
# r1Mod.AICc <- 8 # for AA
cnt_AICc<-sum(apply(delta.AICc[rnks.AICc[,r1Mod.AICc]==1,-r1Mod.AICc] < cutoff.AIC, 1, sum)>0)
cnt_AICc
round(cnt_AICc/Cnt_AICc[1,r1Mod.AICc]*100,1)

#~~~~~~~~~~~~~~~~~~~~~~
# Normalized-MAD statistics of best-performing model
round(mean(minMADs/mean.Nconsumed)*100,2)
round(range(minMADs/mean.Nconsumed)*100,2)

# How many times is a single model the only best model by MAD as judged by being below the data-detependent performance criterion?
cnt_MAD_single<-sum(apply(MADs < cutoff.MAD, 1, sum)==1)
cnt_MAD_single
round(cnt_MAD_single/nrow(MADs)*100,1)

# How many times is a single model the only best model by MAD as judged by being below the *relative* performance criterion?
cnt_MAD_single<-sum(apply(delta.MAD < cutoff.delta.MAD, 1, sum)==1)
cnt_MAD_single
round(cnt_MAD_single/nrow(delta.MAD)*100,1)
nrow(delta.MAD)-cnt_MAD_single
round((nrow(delta.MAD)-cnt_MAD_single)/nrow(delta.MAD)*100,1)

# When the overall top-performing model is best, how often are models other than the overall top-peforming model within the cutoff.delta.MAD performance criterion?
r1Mod.MAD <- which.max(Cnt_MAD[1,])
cnt_MAD <- apply(delta.MAD[rnks.MAD[,r1Mod.MAD ]==1,] < cutoff.delta.MAD, 2, sum)
cnt_MAD
round(100*cnt_MAD/max(cnt_MAD),1)

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
fFirst.AICc<-fSecnd.AICc<-fFirst.MAD<-fSecnd.MAD<-dim(0)
pfFirst.AICc<-pfSecnd.AICc<-pfFirst.MAD<-pfSecnd.MAD<-dim(0)
for(SScut in SScuts){
  
  # Cnt_AICc<-apply(rnks.AICc[sample.sizes>=SScut,],2, # exclude smaller than
  #                 function(x){
  #                   table(factor(x,levels=1:mods.n))})
  Cnt_AICc<-apply(rnks.AICc[sample.sizes<=SScut,],2, # exclude greater than
                  function(x){
                    table(factor(x,levels=1:mods.n))})
  Cnt_AICc <- Cnt_AICc[,mod.order]
  pCnt_AICc <- prop.table(Cnt_AICc,2)
  
  # Cnt_MAD<-apply(rnks.MAD[sample.sizes>=SScut,],2, # exclude smaller than
  #                function(x){
  #                  table(factor(x,levels=1:mods.n))})
  Cnt_MAD<-apply(rnks.MAD[sample.sizes<=SScut,],2, # exclude greater than
                  function(x){
                    table(factor(x,levels=1:mods.n))})
  Cnt_MAD <- Cnt_MAD[,mod.order]
  pCnt_MAD <- prop.table(Cnt_MAD,2)

  fFirst.AICc <- rbind(fFirst.AICc,Cnt_AICc[1,])
  fSecnd.AICc <- rbind(fSecnd.AICc,Cnt_AICc[2,])
  fFirst.MAD <- rbind(fFirst.MAD,Cnt_MAD[1,])
  fSecnd.MAD <- rbind(fSecnd.MAD,Cnt_MAD[2,])
  
  pfFirst.AICc <- rbind(pfFirst.AICc,pCnt_AICc[1,])
  pfSecnd.AICc <- rbind(pfSecnd.AICc,pCnt_AICc[2,])
  pfFirst.MAD <- rbind(pfFirst.MAD,pCnt_MAD[1,])
  pfSecnd.MAD <- rbind(pfSecnd.MAD,pCnt_MAD[2,])

}

# How many datasets are included at the mxSS=300 limit of the plots?
sum(sample.sizes>=mxSS)

# What are stats at the median sample size?
s=median(sample.sizes)
length(sample.sizes[sample.sizes<=s])

fFirst.AICc[which(SScuts==s),]
round(pfFirst.AICc[which(SScuts==s),],3)*100
fSecnd.AICc[which(SScuts==s),]
round(pfSecnd.AICc[which(SScuts==s),],3)*100

fFirst.MAD[which(SScuts==s),]
round(pfFirst.MAD[which(SScuts==s),],3)*100
fSecnd.MAD[which(SScuts==s),]
round(pfSecnd.MAD[which(SScuts==s),],3)*100

# When comparing studies with sample-sizes above versus below the median sample size...
# How many times is a single model the only best model by AICc or MAD as judged by being below the *relative* performance criterion?
# AICc above median:
length(sample.sizes[sample.sizes>s])
cnt_AICc_single<-sum(apply(delta.AICc[sample.sizes>s,] <= cutoff.AIC, 1, sum)==1)
cnt_AICc_single
round(cnt_AICc_single/nrow(delta.AICc[sample.sizes>s,])*100,1)

# AICc at or below median:
length(sample.sizes[sample.sizes<=s])
cnt_AICc_single<-sum(apply(delta.AICc[sample.sizes<=s,] <= cutoff.AIC, 1, sum)==1)
cnt_AICc_single
round(cnt_AICc_single/nrow(delta.AICc[sample.sizes<=s,])*100,1)

# MAD above median:
length(sample.sizes[sample.sizes>s])
cnt_MAD_single<-sum(apply(delta.MAD[sample.sizes>s,] < cutoff.delta.MAD[sample.sizes>s], 1, sum)==1)
cnt_MAD_single
round(cnt_MAD_single/nrow(delta.MAD[sample.sizes>s,])*100,1)

# MAD at or below median:
length(sample.sizes[sample.sizes<=s])
cnt_MAD_single<-sum(apply(delta.MAD[sample.sizes<=s,] < cutoff.delta.MAD[sample.sizes<=s], 1, sum)==1)
cnt_MAD_single
round(cnt_MAD_single/nrow(delta.MAD[sample.sizes<=s,])*100,1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure of top two rankings as a function of sample size
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Darken colors of linear models relative to above
Mcols[c(1,5)]<-'grey40'
ltys <- c(1,1,1,6,2,2,1,1)

pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AICc_and_MAD_toprankBySS.pdf',
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

ylims <- c(0,0.58)

mods<-c(1,5)
legend(xadj1,yadj,legend=colnames(pfFirst.AICc)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       cex=lcex, ncol=1, pt.lwd=0.8, 
       title=expression(italic(k)==1), bty='n')

mods<-c(2,6,7)
legend(xadj2,yadj,legend=colnames(pfFirst.AICc)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==2), bty='n')

mods<-c(3,4,8)
legend(xadj3,yadj,legend=colnames(pfFirst.AICc)[mods],
       lty=ltys[mods], col=Mcols[mods],lwd=1.5,
       pt.cex=1.1,cex=lcex, ncol=1, pt.lwd=0.8,
       title=expression(italic(k)==3), bty='n')

rect(xadj2-0.2,yadj+0.4,xadj2+0.4,yadj-1.5)
text(xadj2+0.02, yadj+0.18, 'Models', adj=0,cex=1)

par(mar=c(2.5,3,0.5,0))
matplot(SScuts,pfFirst.AICc,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='Fraction in first rank',
        xlab='')
  rug(sample.sizes,side=3)
  mtext('A',side=3,line=-1.5,at=mnSS,cex=0.9)
  mtext('AICc',side=3,line=0.3,adj=0.5,cex=1)
  
par(mar=c(3,3,0,0))
matplot(SScuts,pfSecnd.AICc,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='Fraction in second rank',
        # xlab='Sample size greater than...')
        xlab='Sample size less than...')
  rug(sample.sizes,side=3)
  mtext('B',side=3,line=-1.5,at=mnSS,cex=0.9)
  
par(mar=c(2.5,2,0.5,1))
matplot(SScuts,pfFirst.MAD,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='',
        xlab='')
  rug(sample.sizes,side=3)
  mtext('C',side=3,line=-1.5,at=mnSS,cex=0.9)
  mtext('MAD',side=3,line=0.3,adj=0.5,cex=1)
  
par(mar=c(3,2,0,1))
matplot(SScuts,pfSecnd.MAD,las=1,type='l', 
        col=Mcols,lty=ltys,cex=0.5,
        ylim=ylims,lwd=1.5,log='x',
        ylab='',
        # xlab='Sample size greater than...')
        xlab='Sample size less than...')
  rug(sample.sizes,side=3)
  mtext('D',side=3,line=-1.5,at=mnSS,cex=0.9)

dev.off()

