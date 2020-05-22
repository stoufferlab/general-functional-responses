
library(bbmle)
library(RColorBrewer)

# read in the different fits
ffr.fits <- readRDS(
	file='../../../../results/R/OnePredTwoPrey_ffr.fits.Rdata'
)

# re order so the plot is alphabetical top to bottom
ffr.fits <- ffr.fits[rev(1:length(ffr.fits))]

# scrape out the AIC values for the different models
AICs <- t(sapply(
	seq(1,length(ffr.fits)),
	function(x,ffr.fits) {
		unlist(lapply(ffr.fits[[x]]$AICs, function(x){mean(unlist(x))}))
	},
	ffr.fits=ffr.fits
))
colnames(AICs) <- c(
	"H1",
	"H2.SS",
	"H2.SG",
	"H2.GS",
	"H2.GG",
	"H2.HHI",
	"H2.HHE"
)
AICs <- as.data.frame(AICs)

# break a tie between the two hybrid-hybrid models
AICs[,"H2.HH"] <- pmin(AICs[,"H2.HHI"],AICs[,"H2.HHE"])
AICs[,"H2.HHI"] <- AICs[,"H2.HHE"] <- NULL

# DEBUG
# remove the SG and GS models?
AICs <- AICs[,c("H1","H2.SS","H2.GG","H2.HH")]

labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
labels <- paste(labels,"  ")

# determine things based on deltaAIC
minAICs <- apply(AICs, 1, min)
dAICs <- AICs - minAICs
dAICs[dAICs<2] <- 0
rnkAICs <- t(apply(dAICs, 1, rank, ties.method='first'))
colnames(rnkAICs) <- colnames(AICs)

# Define delta AICc cut-off for "indistinguishably well performing" models
delAICcutoff <- 2

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~

colnames(AICs)
CR<-brewer.pal(n = 7, name = 'Blues')
Mcols <- CR[c(2,4,6)]
# Mcols <- c("white",Mcols)
Mcols <- c(Mcols, brewer.pal(n = 9, name = 'Reds')[7])
Mpch <- c(rep(21,7),rep(22,4))

# generate the figure

pdf(
    '../../../../results/R/OnePredTwoPrey_figs/OnePredTwoPrey_AIC_ranks.pdf',
    height=2.875,
    width=2.25
)

par(
    mar=c(3,7,2.5,0.5),
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
    axis(2, at=1:nrow(rnkAICs), labels=labels, cex.axis=0.5, las=2, lwd=0, lwd.ticks=1)
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)

    # Which models have delta-AIC within X=2 of best-performing model?
    xats <-table(which(dAICs < delAICcutoff, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))

    # segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    # segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
  
    # shade behind ties
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.4
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
        -0.5,nrow(rnkAICs)+6,legend=colnames(rnkAICs),
        pch=Mpch, pt.bg=Mcols, col='black', bg='white',
        horiz=FALSE, pt.cex=1, cex=0.6, ncol=4, title='Model',
        # xjust=0.5,
        title.adj=0.29,
        bty='n'
    )
    par(xpd=FALSE)
dev.off()


# # profile all of the datasets
# ffr.cfs <- lapply(
# 	ffr.fits,
# 	function(x,modeltype,parameters){
# 		if(x$estimates[[modeltype]]["n",1,1] != 1){
# 			return(NULL)
# 		}else{
# 			foobar <- list()
# 			foobar[[modeltype]] <- try(
# 				confint(
# 					proffun(
# 						x$fits[[modeltype]],
# 						which=parameters,
# 						try_harder=TRUE,
# 						level=0.68,
# 						tol.newmin=Inf
# 					)
# 				)
# 			)
# 			foobar
# 		}
# 	},
# 	modeltype=modeltype,
# 	parameters=parameters
# )

# # save the mega container which includes all FR fits
# save(ffr.cfs,file='../../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata')
