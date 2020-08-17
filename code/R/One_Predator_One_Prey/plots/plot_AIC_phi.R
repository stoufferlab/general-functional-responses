source('../../lib/profile_coefs.R')
source('../../lib/plot_coefs.R')
source('../../lib/depletion_check.R') 
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

library(RColorBrewer)
library(Hmisc) # for LaTeX table export
library(stringr)
options(xdvicmd='open')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Stouffer.Novak.I", order.parm="phi_denom")
ffr.fits <- ffr.fits[fit.order]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))

# need to doctor a label to match the supps
labels[which(labels=="Long_2012")] <- "Long_2012b"

# replace underscores
labels<-gsub('_',' ',labels)

# Grab summary of AIC estimates across bootstrapped fits
stat <- '50%' # use "50%" or "mean"

AIC.H1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.I'][[1]][stat]}))
AIC.H2 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Holling.II'][[1]][stat]}))
AIC.BD <- unlist(lapply(ffr.fits, function(x){ x$AIC['Beddington.DeAngelis'][[1]][stat]}))
AIC.CM <- unlist(lapply(ffr.fits, function(x){ x$AIC['Crowley.Martin'][[1]][stat]}))
AIC.SN1 <- unlist(lapply(ffr.fits, function(x){ x$AIC['Stouffer.Novak.I'][[1]][stat]}))

AICs <- data.frame(AIC.H1, AIC.H2, AIC.BD, AIC.CM, AIC.SN1)
colnames(AICs) <- sub('AIC.', '', colnames(AICs))
colnames(AICs)[5] <- "G"
rownames(AICs) <- labels

# Define color of each model
colnames(AICs)
CR<-brewer.pal(n = 9, name = 'Blues')
Mcols <- c('white',CR[c(3,6,9)])
Mcols <- c(Mcols, brewer.pal(n = 9, name = 'Reds')[7])
Mpch <- c(rep(21,4),rep(22,4))

minAICs <- apply(AICs, 1, min)
dAICs <- AICs - minAICs
# dAICs[dAICs<2] <- 0

rnkAICs <- t(apply(dAICs, 1, rank, ties.method='first'))
colnames(rnkAICs) <- colnames(AICs)

# Define delta AIC cut-off for "indistinguishably well performing" models
delAICcutoff <- 2

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_AIC_phi.pdf',
    height=6,
    width=4.5
)

layout(mat = matrix(
        c(1, 2), 
        nrow = 1, 
        ncol = 2
       ),
       heights = c(1),
       widths = c(3.25, 3.5)
)

par(
    mar=c(2.5,2.75,0.5,5.75),
    oma = c(0, 0, 0, 0),
    mgp=c(1.25,0.1,0.0),
    tcl=-0.1,
    las=1,
    cex=0.7,
    yaxs='i'
)

    plot(1:nrow(rnkAICs), 1:nrow(rnkAICs),
         type='n', yaxt='n',
         xlim=c(0.5,ncol(rnkAICs)+0.5),
         ylim=c(0,nrow(rnkAICs)+1),
         xlab='Model rank by AIC',
         ylab='',
         axes=F)
    # rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30

# add some padding to pretty everything up
labels <- str_pad(labels, side="both", width=max(str_length(labels))+2)

    axis(
        4,
        at=1:nrow(rnkAICs),
        labels=labels,
        cex.axis=0.5,
        las=2,
        lwd=0,
        lwd.ticks=1,
        hadj=0.5,
        mgp=c(0,3,0)
    )
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)

    # Which models have delta-AIC within X=2 of best-performing model?
    xats <-table(which(dAICs < delAICcutoff, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))

    # segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    # segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
  
    # shade behind ties
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.50

    pyats[1:2] <- pyats[1:2] - 0.5
    pyats[(length(pyats)-1):length(pyats)] <- pyats[(length(pyats)-1):length(pyats)]+0.5

    polygon(pxats,pyats,col=grey(0.666),border=NA)
    
    for(m in 1:ncol(rnkAICs)){
      points(rnkAICs[,m], 1:nrow(rnkAICs), 
             type='p',  col='black', 
             bg=Mcols[m], pch=Mpch[m],
             cex=1, lwd=0.2)
    }  
    box(lwd=1)
    par(xpd=TRUE)
    legend(
        -0.75,nrow(rnkAICs)/2,legend=colnames(rnkAICs),
        xjust=0.5,
        pch=Mpch, pt.bg=Mcols, col='black', bg='white',
        horiz=FALSE, pt.cex=1,cex=0.6, ncol=1, title='Model',
        yjust = 0.5,
        bty='n'
    )
    par(xpd=FALSE)

# fits for which the SN1 model is not supported are colored red
color.vector <- numeric(length(ffr.fits))
for(i in 1:length(ffr.fits)){
  if(dAICs$G[i] > 2){
    color.vector[i] <- CR[7]
  }else{
    if(rowSums(dAICs[i,]<2)==1){
      color.vector[i] <- "black"
    }else{
      color.vector[i] <- "black"
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ffr.fits <- profile_coefs(
  ffr.fits,
  model='Stouffer.Novak.I',
  point.est='median',
  printWarnings = FALSE,
  which.pars = "phi_denom"
)
save(ffr.fits, file='../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.SN1.Rdata')
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels <- rep("",length(labels))

par(mar=c(2.5,0.25,0.5,0.5))

plot.coefs(
    ffr.fits,
    model="Stouffer.Novak.I",
    parameter="phi_denom",
    ilink=identity,
    point.est='median',
    plot.SEs=TRUE,
    display.outlier.ests=FALSE,
    xlab="",
    labels=labels,
    vertLines = c(0,1),
    xlim=c(-3,3),
    ylim=c(0,length(ffr.fits)+1),
    color.vector = color.vector,
    SE.lty=c(1,6,3)
)

mtext(
  expression(plain('Effect of feeding on interference, ')*italic(varphi[I[k]][F[ki]])),
  1,
  cex=0.7,
  line=1.5
)

dev.off()
