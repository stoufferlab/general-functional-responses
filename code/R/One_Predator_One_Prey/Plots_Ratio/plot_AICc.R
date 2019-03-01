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

# Grab summary of AICc estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

AICc.H1 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.Type.I'][[1]][stat]}))
AICc.H2 <- unlist(lapply(ffr.fits, function(x){ x$AICc['Holling.Type.II'][[1]][stat]}))
AICc.R <- unlist(lapply(ffr.fits, function(x){ x$AICc['Ratio'][[1]][stat]}))
AICc.AG <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Ginzburg'][[1]][stat]}))
AICc.HV <- unlist(lapply(ffr.fits, function(x){ x$AICc['Hassell.Varley'][[1]][stat]}))
AICc.AA <- unlist(lapply(ffr.fits, function(x){ x$AICc['Arditi.Akcakaya'][[1]][stat]}))
AICc.BD <- unlist(lapply(ffr.fits, function(x){ x$AICc['Beddington.DeAngelis'][[1]][stat]}))
AICc.CM <- unlist(lapply(ffr.fits, function(x){ x$AICc['Crowley.Martin'][[1]][stat]}))


AICcs <- data.frame(AICc.H1, AICc.H2, AICc.R, AICc.AG, AICc.HV, AICc.AA, AICc.BD, AICc.CM)

rownames(AICcs)<-sub('..Dataset_Code.','', rownames(AICcs))
rownames(AICcs)<-sub('.R.50%','', rownames(AICcs))
rownames(AICcs)<-sub('.R.mean','', rownames(AICcs))
rownames(AICcs) <- paste0(rownames(AICcs), ' (',SS,')')
colnames(AICcs) <- sub('AICc.', '', colnames(AICcs))

minAICcs <- apply(AICcs, 1, min)

dAICcs <- AICcs - minAICcs

rnkAICcs <- t(apply(dAICcs, 1, rank, ties.method='first'))
colnames(rnkAICcs) <- colnames(AICcs)

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
Mcols <- brewer.pal(ncol(rnkAICcs),'Set3')

pdf('../../../../results/R/OnePredOnePrey_AICc_ranks.pdf',height=6,width=2.25)
par(mar=c(3,9,2.5,0.5), mgp=c(1.5,0.2,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
  plot(1:nrow(rnkAICcs), 1:nrow(rnkAICcs),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnkAICcs)),
       ylim=c(0,nrow(rnkAICcs)+1),
       xlab='Model rank by AICc',
       ylab='',
       axes=F)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey30")
  axis(2, at=1:nrow(rnkAICcs), labels=rownames(rnkAICcs), cex.axis=0.5, las=2)
  axis(1, cex.axis=0.7, mgp=c(1.25,0,0))
  box(lwd=1)
  

  # Which models have delta-AICc within X of best-performing model?
  # (Could use either of the next sections depending on preference)
  xats <-table(which(dAICcs<2,arr.ind=T)[,1])+0.5

  xats <-table(which(dAICcs<3,arr.ind=T)[,1])+0.5
  yats <- 0:(length(xats))+0.5
  segments(xats,yats[-length(yats)],xats,yats[-1],col='white')
  segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='white')

  pxats<-c(1,rep(xats,each=2),1)
  pyats<-rep(0:(length(xats)),each=2)+0.5
  polygon(pxats,pyats,col='grey40',border=NA)
  
  for(m in 1:ncol(rnkAICcs)){
    points(rnkAICcs[,m], 1:nrow(rnkAICcs), type='p', pch=21, col=Mcols[m], bg=Mcols[m])
  }  
  par(xpd=TRUE)
  legend(-8,nrow(rnkAICcs)+6,legend=colnames(rnkAICcs),pch=21, pt.bg=Mcols, col='black', bg='white', horiz=T,pt.cex=1.1,cex=0.6, ncol=2, title='Model')
dev.off()

#~~~~~
# dAICc
#~~~~~
plot(1:nrow(dAICcs), 1:nrow(dAICcs),
     type='n', yaxt='n',
     # xlim=c(0,max(dAICcs)),
     xlim=c(0,4),
     xlab='Model delta-AICc',
     ylab=''
)
axis(side=2, at=1:nrow(dAICcs), labels=rownames(dAICcs), cex.axis=0.5, las=2)

Mcols <- rainbow(ncol(dAICcs))
for(m in 1:ncol(dAICcs)){
  points(dAICcs[,m], 1:nrow(dAICcs), pch=21, col=Mcols[m], bg=Mcols[m])
}  

legend('bottomright',legend=colnames(dAICcs),pch=21, pt.bg=Mcols, col=Mcols)

