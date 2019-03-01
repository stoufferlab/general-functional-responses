source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

load('../../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

library(RColorBrewer)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
SS <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))


# Grab summary of AIC estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

AIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.Type.I'][[1]][stat]}))
AIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.Type.II'][[1]][stat]}))
AIC.R <- unlist(lapply(ffr.fits, function(x){ x$AIC['Ratio'][[1]][stat]}))
AIC.AG <- unlist(lapply(ffr.fits, function(x){ x$AIC['Arditi.Ginzburg'][[1]][stat]}))
AIC.HV <- unlist(lapply(ffr.fits, function(x){ x$AIC['Hassell.Varley'][[1]][stat]}))
AIC.AA <- unlist(lapply(ffr.fits, function(x){ x$AIC['Arditi.Akcakaya'][[1]][stat]}))
AIC.BD <- unlist(lapply(ffr.fits, function(x){ x$AIC['Beddington.DeAngelis'][[1]][stat]}))
AIC.CM <- unlist(lapply(ffr.fits, function(x){ x$AIC['Crowley.Martin'][[1]][stat]}))


AICs <- data.frame(AIC.H1, AIC.H2, AIC.R, AIC.AG, AIC.HV, AIC.AA, AIC.BD, AIC.CM)

rownames(AICs)<-sub('..Dataset_Code.','', rownames(AICs))
rownames(AICs)<-sub('.R.50%','', rownames(AICs))
rownames(AICs)<-sub('.R.mean','', rownames(AICs))
rownames(AICs) <- paste0(rownames(AICs), ' (',SS,')')
colnames(AICs) <- sub('AIC.', '', colnames(AICs))

minAICs <- apply(AICs, 1, min)

dAICs <- AICs - minAICs

rnkAICs <- t(apply(dAICs, 1, rank, ties.method='first'))
colnames(rnkAICs) <- colnames(AICs)

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
Mcols <- brewer.pal(ncol(rnkAICs),'Set3')

pdf('../../../../results/R/OnePredOnePrey_AIC_ranks.pdf',height=6,width=2.25)
par(mar=c(3,9,2.5,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
  plot(1:nrow(rnkAICs), 1:nrow(rnkAICs),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnkAICs)),
       ylim=c(0,nrow(rnkAICs)+1),
       xlab='Model rank by AIC',
       ylab='',
       axes=F)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey30")
  axis(2, at=1:nrow(rnkAICs), labels=rownames(rnkAICs), cex.axis=0.5, las=2)
  axis(1, cex.axis=0.7, mgp=c(1.25,0,0))
  box(lwd=1)
  
  xats <-table(which(dAICs<3,arr.ind=T)[,1])+0.5
  yats <- 0:(length(xats))+0.5
  segments(xats,yats[-length(yats)],xats,yats[-1],col='white')
  segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='white')

  pxats<-c(1,rep(xats,each=2),1)
  pyats<-rep(0:(length(xats)),each=2)+0.5
  polygon(pxats,pyats,col='grey40',border=NA)
  
  for(m in 1:ncol(rnkAICs)){
    points(rnkAICs[,m], 1:nrow(rnkAICs), type='p', pch=21, col=Mcols[m], bg=Mcols[m])
  }  
  par(xpd=TRUE)
  legend(-8,nrow(rnkAICs)+6,legend=colnames(rnkAICs),pch=21, pt.bg=Mcols, col='black', bg='white', horiz=T,pt.cex=1.1,cex=0.6, ncol=2, title='Model')
dev.off()

#~~~~~
# dAIC
#~~~~~
plot(1:nrow(dAICs), 1:nrow(dAICs),
     type='n', yaxt='n',
     # xlim=c(0,max(dAICs)),
     xlim=c(0,4),
     xlab='Model delta-AIC',
     ylab=''
)
axis(side=2, at=1:nrow(dAICs), labels=rownames(dAICs), cex.axis=0.5, las=2)

Mcols <- rainbow(ncol(dAICs))
for(m in 1:ncol(dAICs)){
  points(dAICs[,m], 1:nrow(dAICs), pch=21, col=Mcols[m], bg=Mcols[m])
}  

legend('bottomright',legend=colnames(dAICs),pch=21, pt.bg=Mcols, col=Mcols)

