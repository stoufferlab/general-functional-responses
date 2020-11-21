
# libraries that may or may not be needed below
library(RColorBrewer)
library(shape)
library(stringr)

# generate the AIC table
source('../summaries/make_AICc_table.R')

# generate the BIC table
source('../summaries/make_BIC_table.R')

# read in the profiled CIs for the phi values
load(
  file='../../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata'
)

# make sure things are named correctly
names(ffr.cfs) <- unlist(lapply(ffr.cfs, function(x) x$study.info$datasetName))

# read in the cheat sheet to figure out which HH model between the two parameter transformations
source('which_model.R')

# figure out the average phi value across datasets
mean_phi <- sapply(
  rownames(rnkAICcs),
  function(x,ffr.cfs,which.model){
    mean(ffr.cfs[[x]]$profiles[[paste0("Holling II Hybrid Hybrid ",which.model[x])]]$est)
  },
  ffr.cfs=ffr.cfs,
  which.model=which.model
)

# reorder the AICc and BIC tables and the profiled fits
dAICcs <- dAICcs[order(mean_phi),]
rnkAICcs <- rnkAICcs[order(mean_phi),]
dBICs <- dBICs[order(mean_phi),]
rnkBICs <- rnkBICs[order(mean_phi),]
ffr.cfs <- ffr.cfs[order(mean_phi)]

# colors and symbols for different models
CR<-brewer.pal(n = 9, name = 'Blues')
Mcols <- CR[c(6,8)]
Mcols <- c("white",Mcols)
Mcols <- c(Mcols, brewer.pal(n = 9, name = 'Reds')[7])
Mpch <- c(rep(21,3),rep(22,4))

# doctor one label so that it matches the table in the supps
labels <- rownames(rnkAICcs)
labels[which(labels=="Long_2012_2p")] <- "Long_2012a"

# create the dataset labels for the figures
labels<-gsub('_',' ',labels)
labels <- str_pad(labels, side="both", width=max(str_length(labels))+2)

# generate the figure
filename <- '../../../../results/R/OnePredTwoPrey_figs/OnePredTwoPrey_AICc_and_BIC_ranks'
# pdf(
#     paste0(filename,'.pdf'),
#     height=2.875,
#     width=4.25
# )
setEPS()
postscript(paste0(filename,'.eps'),
           height=2.875,
           width=4.25
)

layout(mat = matrix(
        c(1, 2), 
        nrow = 1, 
        ncol = 2
       ),
       heights = c(1),
       widths = c(2.75,2.75)
)

par(
    mar=c(2.5,2.75,0.5,4.65),
    oma = c(0, 2.5, 0, 0),
    mgp=c(1.25,0.1,0.0),
    tcl=-0.1,
    las=1,
    cex=0.7,
    yaxs='i'
)

#########################################################
# AICc plot
#########################################################

plot(1:nrow(rnkAICcs), 1:nrow(rnkAICcs),
     type='n', yaxt='n',
     xlim=c(0.5,ncol(rnkAICcs)+0.5),
     ylim=c(0,nrow(rnkAICcs)+1),
     xlab='Model rank by AICc',
     ylab='',
     axes=F
)
    # rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
axis(
  4,
  at=1:nrow(rnkAICcs),
  labels=labels,
  cex.axis=0.5,
  las=2,
  lwd=0,
  lwd.ticks=1,
  hadj=0.5,
  mgp=c(0,2.5,0)
)
axis(1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)

# Which models have delta-AICc within X=2 of best-performing model?
xats <-table(which(dAICcs < delAICccutoff, arr.ind=T)[,1])+0.5
yats <- 0:(length(xats))

# segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
# segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')

# shade behind ties
pxats<-c(0,rep(xats,each=2),0)
pyats<-rep(0:(length(xats)),each=2)+0.5
pyats[1:2] <- pyats[1:2]-0.5
pyats[(length(pyats)-1):length(pyats)] <- pyats[(length(pyats)-1):length(pyats)]+0.5
polygon(pxats,pyats,col=grey(0.875),border=NA)

for(m in 1:ncol(rnkAICcs)){
  points(rnkAICcs[,m], 1:nrow(rnkAICcs), 
         type='p',  col='black', 
         bg=Mcols[m], pch=Mpch[m],
         cex=1, lwd=0.2
  )
}

box(lwd=1)
par(xpd=TRUE)
legend(
    -0.75,nrow(rnkAICcs)/2,
    legend=c(
      'H1',
      "H2",
      "H2-m",
      "H2-g" #expression(paste("H2-",phi))
    ),
    pch=Mpch, pt.bg=Mcols, col='black', bg='white',
    horiz=FALSE, pt.cex=1, cex=0.6, ncol=1, title='  Model',
    xjust=0.5, yjust=0.5,
    title.adj=0.29,
    bty='n'
)
par(xpd=FALSE)

#########################################################
# BIC plot
#########################################################

par(
  mar=c(2.5,0.25,0.5,7.1250)
  # oma = c(0, 0, 0, 2)
)

plot(1:nrow(rnkBICs), 1:nrow(rnkBICs),
     type='n', yaxt='n',
     xlim=c(0.5,ncol(rnkBICs)+0.5),
     ylim=c(0,nrow(rnkBICs)+1),
     xlab='Model rank by BIC',
     ylab='',
     axes=F
)
    # rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = "white") # grey30
axis(
  2,
  at=1:nrow(rnkBICs),
  labels=rep("",nrow(rnkBICs)), #labels,
  cex.axis=0.5,
  las=2,
  lwd=0,
  lwd.ticks=1,
  hadj=0.5,
  mgp=c(0,2.5,0)
)
axis(1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)

# Which models have delta-BIC within X=2 of best-performing model?
xats <-table(which(dBICs < delBICcutoff, arr.ind=T)[,1])+0.5
yats <- 0:(length(xats))

# segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
# segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')

# shade behind ties
pxats<-c(0,rep(xats,each=2),0)
pyats<-rep(0:(length(xats)),each=2)+0.5
pyats[1:2] <- pyats[1:2]-0.5
pyats[(length(pyats)-1):length(pyats)] <- pyats[(length(pyats)-1):length(pyats)]+0.5
polygon(pxats,pyats,col=grey(0.875),border=NA)

for(m in 1:ncol(rnkBICs)){
  points(rnkBICs[,m], 1:nrow(rnkBICs), 
         type='p',  col='black', 
         bg=Mcols[m], pch=Mpch[m],
         cex=1, lwd=0.2
  )
}

box(lwd=1)

dev.off()
