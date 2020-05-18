
library(bbmle)
library(RColorBrewer)
library(stringr)

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
labels <- str_pad(labels, side="both", width=max(str_length(labels))+2)

# labels <- paste(labels,"  ")

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
CR<-brewer.pal(n = 9, name = 'Blues')
Mcols <- CR[c(6,8)]
Mcols <- c("white",Mcols)
Mcols <- c(Mcols, brewer.pal(n = 9, name = 'Reds')[7])
Mpch <- c(rep(21,7),rep(22,4))

# generate the figure

pdf(
    '../../../../results/R/OnePredTwoPrey_figs/OnePredTwoPrey_AIC_phi.pdf',
    height=2.875,
    width=3.25
)

layout(mat = matrix(
        c(1, 2), 
        nrow = 1, 
        ncol = 2
       ),
       heights = c(1),
       widths = c(2.5, 2.0)
)

par(
    mar=c(2.5,2.75,0.5,4.65),
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
    axis(
      4,
      at=1:nrow(rnkAICs),
      labels=labels,
      cex.axis=0.5,
      las=2,
      lwd=0,
      lwd.ticks=1,
      hadj=0.5,
      mgp=c(0,2.5,0)
    )
    axis(1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)

    # Which models have delta-AIC within X=2 of best-performing model?
    xats <-table(which(dAICs < delAICcutoff, arr.ind=T)[,1])+0.5
    yats <- 0:(length(xats))

    # segments(xats,yats[-length(yats)],xats,yats[-1],col='black')
    # segments(xats[-length(xats)],yats+1,xats[-1],yats+1,col='black')
  
    # shade behind ties
    pxats<-c(0,rep(xats,each=2),0)
    pyats<-rep(0:(length(xats)),each=2)+0.5
    pyats[1:2] <- pyats[1:2]-0.5
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
        -0.75,nrow(rnkAICs)/2,
        legend=c(
          'H1',
          "H2 (S)",
          "H2 (G)",
          expression(paste("H2 (",phi,")"))
        ),
        pch=Mpch, pt.bg=Mcols, col='black', bg='white',
        horiz=FALSE, pt.cex=1, cex=0.6, ncol=1, title='  Model',
        xjust=0.5, yjust=0.5,
        title.adj=0.29,
        bty='n'
    )
    par(xpd=FALSE)

#########################################################

# read in the profiled CIs
load(
  file='../../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata'
)

# reverse order
ffr.cfs <- ffr.cfs[rev(1:length(ffr.cfs))]

# x axis limits
xlim <- c(-3,4.5)

# some colors
CR<-brewer.pal(n = 9, name = 'Blues')

# jitter
dy <- 0.4

# dataset labels
labels <- unlist(lapply(ffr.cfs, function(x) x$study.info$datasetName))

# cheatsheet for which model to plot
# NOTE: this is done by hand and shameful
which.model <- data.frame(rbind(
  c("Chan_2017_c",1),
  c("Chan_2017_l",2),
  c("Colton_1987_1st24",1),
  c("Colton_1987_2nd24",1),
  c("Elliot_2006_Instar2",1),
  c("Elliot_2006_Instar3",2),
  c("Elliot_2006_Instar4",1),
  c("Elliot_2006_Instar4Baet",1),
  c("Elliot_2006_Instar5",1),
  c("Elliot_2006_Instar5Baet",1),
  c("Iyer_1996_Bc",1),
  c("Iyer_1996_Bp",1),
  c("Iyer_1996_Br",1),
  c("Kalinkat_2011_Anch",1),
  c("Kalinkat_2011_Cal",2),
  c("Kalinkat_2011_Harp",1),
  c("Kalinkat_2011_Pard",2),
  c("Kalinkat_2011_Troch",2),
  c("Krylov_1992ii",1),
  c("Lester_2002_Af_duets",1),
  c("Lester_2002_Af_eggs",1),
  c("Lester_2002_Ty_duets",1),
  c("Lester_2002_Ty_eggs",1),
  c("Long_2012b",1),
  c("Nachappa_2006",2),
  c("Ranta_1985_10",1),
  c("Ranta_1985_13",1),
  c("Ranta_1985_18",1),
  c("Ranta_1985_Ad",2),
  c("Wong_2005_rc",1),
  c("Wong_2005_ss",2)
))
rownames(which.model) <- which.model[,1]
which.model[,1] <- NULL

par(mar=c(2.5,0.25,0.5,0.5))

# plot.coefs <- function(ffr.fits, modeltype, parameters, plot.SEs=FALSE, ffr.cfs=NULL, ilink=identity, xlim=NULL, ... ){
  # scaffold of a plot that doesn't actually show anything
  plot(
    y = 1:31,
    x = 1:31,
    type='n',
    yaxt='n',
    xlim=xlim,
        ylim=c(0,32),
        xlab='Effect of feeding on feeding',
        ylab='',
        axes=FALSE
  )

  # mark where the existing models fall
  abline(v=c(0,1), lty=2, col="grey")
  # abline(v=1, lty=2)

  # tidy up the y axis labels
  labels<-gsub('_',' ',labels)
  labels <- paste(labels,"  ")

  # tick marks to indicate different data sets
  axis(side=2, at=1:length(ffr.cfs), labels=rep("",length(ffr.cfs)), cex.axis=0.5, las=1, lwd=0, lwd.ticks=1)
  axis(side=1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)
  box(lwd=1)

  # length of arrows to indicate values off the plot
  delta.arrow <- 0.04*diff(xlim)

  # plot these bad boys
  for(i in 1:length(ffr.cfs)){
    x <- ffr.cfs[[i]]
    dname <- x$study.info$datasetName

    col <- ifelse(
      dAICs[i,"H2.HH"]<=2,
      "black",
      CR[7]
    )
    for(j in 1:nrow(x$profiles[[1]])){
      # the median estimate is easy to determine regardless of the type of data
      mm <- x$profiles[[which.model[dname,1]]][j,"est"]
      lb <- x$profiles[[which.model[dname,1]]][j,"lb"]
      ub <- x$profiles[[which.model[dname,1]]][j,"ub"]

      # make all lines the equivalent for now
      lty <- switch(x$profiles[[which.model[dname,1]]][j,"method"],
        profile = 1,
        bootstrap = 3,
        quadratic = 6
      )

      lb <- ifelse(is.finite(lb), lb, xlim[1])
      ub <- ifelse(is.finite(ub), ub, xlim[2])

      lb <- ifelse(lb < xlim[1], xlim[1], lb)
      ub <- ifelse(ub > xlim[2], xlim[2], ub)

      if(mm < xlim[1]){
        Arrows(
              xlim[1]+delta.arrow,
              i+dy*(j-1.5),
              xlim[1]-0.10,
              i+dy*(j-1.5),
              segment=FALSE,
              arr.type="triangle",
              arr.adj=1,
              arr.length=0.05,
              arr.width=0.025,
              arr.col=col,
              lcol=col,
              lty=1
            )
            segments(
          xlim[1]+delta.arrow,
          i+dy*(j-1.5),
          xlim[1]-0.15,
          i+dy*(j-1.5),
          col=col,
          lty=lty,
          lwd=0.5
        )

        text(
          xlim[1]+0.15*delta.arrow,
          i+dy*(j-1.5),
          format.number(mm),
          pos=4,
          cex=0.3
        )
      }else if(mm > xlim[2]){
        Arrows(
              xlim[2]-delta.arrow,
              i+dy*(j-1.5),
              xlim[2]+0.10,
              i+dy*(j-1.5),
              segment=FALSE,
              arr.type="triangle",
              arr.adj=1,
              arr.length=0.05,
              arr.width=0.025,
              arr.col=col,
              lcol=col,
              lty=1
            )
            segments(
          xlim[2]-delta.arrow,
          i+dy*(j-1.5),
          xlim[2]+0.15,
          i+dy*(j-1.5),
          col=col,
          lty=lty,
          lwd=0.5
        )
        text(
          xlim[2]-0.15*delta.arrow,
          i+dy*(j-1.5),
          format.number(mm),
          pos=2,
          cex=0.3
        )
      }else{
        # if(lb > xlim[1] & ub < xlim[2]){
        # draw the error bars

        segments(lb, i+dy*(j-1.5), ub, i+dy*(j-1.5), col=col, lty=lty, lwd=0.5)

        # arrow up the limiting cases
        if(lb==xlim[1]){
          Arrows(
                xlim[1]+delta.arrow,
                i+dy*(j-1.5),
                xlim[1]-0.10,
                i+dy*(j-1.5),
                segment=FALSE,
                arr.type="triangle",
                arr.adj=1,
                arr.length=0.05,
                arr.width=0.025,
                arr.col=col,
                lcol=col,
                lty=1
              )
          # arrows(xlim[1], i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        }
        if(ub==xlim[2]){
          Arrows(
                xlim[2]-delta.arrow,
                i+dy*(j-1.5),
                xlim[2]+0.10,
                i+dy*(j-1.5),
                segment=FALSE,
                arr.type="triangle",
                arr.adj=1,
                arr.length=0.05,
                arr.width=0.025,
                arr.col=col,
                lcol=col,
                lty=1
              )
        }
        
        # plot the actual mean estimate
        points(y=i+dy*(j-1.5),x=mm,col=col,bg=col,pch=21,cex=0.5)
      }
      # }else{
        # if(mm[j] > xlim[2]){
        #   arrows(xlim[2]-delta.arrow, i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        # }else{
        #   arrows(xlim[1]+delta.arrow, i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        # }
      # }
    }
  }

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
