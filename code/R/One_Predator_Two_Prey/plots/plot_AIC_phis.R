
# libraries that may or may not be needed below
library(RColorBrewer)
library(shape)
library(stringr)

# generate the AIC tables
source('../results/make_AIC_table.R')

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
  rownames(rnkAICs),
  function(x,ffr.cfs,which.model){
    mean(ffr.cfs[[x]]$profiles[[paste0("Holling II Hybrid Hybrid ",which.model[x])]]$est)
  },
  ffr.cfs=ffr.cfs,
  which.model=which.model
)

# reorder the AIC tables and the profiled fits
dAICs <- dAICs[order(mean_phi),]
rnkAICs <- rnkAICs[order(mean_phi),]
ffr.cfs <- ffr.cfs[order(mean_phi)]

# colors and symbols for different models
CR<-brewer.pal(n = 9, name = 'Blues')
Mcols <- CR[c(6,8)]
Mcols <- c("white",Mcols)
Mcols <- c(Mcols, brewer.pal(n = 9, name = 'Reds')[7])
Mpch <- c(rep(21,3),rep(22,4))

# doctor one label so that it matches the table in the supps
labels <- rownames(rnkAICs)
labels[which(labels=="Long_2012_2p")] <- "Long_2012a"

# create the dataset labels for the figures
labels<-gsub('_',' ',labels)
labels <- str_pad(labels, side="both", width=max(str_length(labels))+2)

# generate the figure
pdf(
    '../../../../results/R/OnePredTwoPrey_figs/OnePredTwoPrey_AIC_phi.pdf',
    height=2.875,
    width=4.25
)

layout(mat = matrix(
        c(1, 2), 
        nrow = 1, 
        ncol = 2
       ),
       heights = c(1),
       widths = c(2.75,3.25)
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

#########################################################
# AIC plot
#########################################################

plot(1:nrow(rnkAICs), 1:nrow(rnkAICs),
     type='n', yaxt='n',
     xlim=c(0.5,ncol(rnkAICs)+0.5),
     ylim=c(0,nrow(rnkAICs)+1),
     xlab='Model rank by AIC',
     ylab='',
     axes=F
)
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
         cex=1, lwd=0.2
  )
}

box(lwd=1)
par(xpd=TRUE)
legend(
    -0.75,nrow(rnkAICs)/2,
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
# phi plot
#########################################################

# x axis limits
xlim <- c(-3,4.5)

# jitter
dy <- 0.4

# plot margins
par(mar=c(2.5,0.25,0.5,0.5))

# scaffold of a plot that doesn't actually show anything
plot(
  y = 1:length(ffr.cfs),
  x = 1:length(ffr.cfs),
  type='n',
  yaxt='n',
  xlim=xlim,
      ylim=c(0,1+length(ffr.cfs)),
      xlab="",
      ylab='',
      axes=FALSE
)

mtext(
  expression(plain('Effects of feeding on feeding, ')*italic(varphi[F[ki]][F[kj]])*plain(' and ')*italic(varphi[F[kj]][F[ki]])),
  1,
  cex=0.7,
  line=1.6
)

# mark where the existing models fall
abline(v=c(0,1), lty=2, col="grey")
# abline(v=1, lty=2)

# tick marks to indicate different data sets
axis(side=2, at=1:length(ffr.cfs), labels=rep("",length(labels)), cex.axis=0.5, las=1, lwd=0, lwd.ticks=1)
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
      model.name <- paste0("Holling II Hybrid Hybrid ",which.model[dname])

      # the median estimate is easy to determine regardless of the type of data
      mm <- x$profiles[[model.name]][j,"est"]
      lb <- x$profiles[[model.name]][j,"lb"]
      ub <- x$profiles[[model.name]][j,"ub"]

      # make all lines the equivalent for now
      lty <- switch(x$profiles[[model.name]][j,"method"],
        profile = 1,
        bootstrap = 3,
        quadratic = 6
      )

      # if we have NA bounds then just make them "maximal"
      lb <- ifelse(is.finite(lb), lb, xlim[1])
      ub <- ifelse(is.finite(ub), ub, xlim[2])

      # lb <- ifelse(lb < xlim[1], xlim[1], lb)
      # ub <- ifelse(ub > xlim[2], xlim[2], ub)

      if(ub < xlim[1]){
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
      }else if(lb > xlim[2]){
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

        segments(
          max(lb,xlim[1]),
          i+dy*(j-1.5),
          min(ub,xlim[2]),
          i+dy*(j-1.5),
          col=col,
          lty=lty,
          lwd=0.5
        )

        # arrow up the limiting cases
        if(lb<=xlim[1]){
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
        if(ub>=xlim[2]){
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
        
        par(xpd=TRUE)
        # plot the actual mean estimate
        if(mm > xlim[1] && mm < xlim[2]){
          points(y=i+dy*(j-1.5),x=mm,col=col,bg=col,pch=21,cex=0.5)
        }else if(mm > xlim[2]){
          points(y=i+dy*(j-1.5),x=(0*xlim[2]+2*par("usr")[2])/2,col=col,bg=col,pch=21,cex=0.5)
        }else{
          points(y=i+dy*(j-1.5),x=(0*xlim[1]+2*par("usr")[1])/2,col=col,bg=col,pch=21,cex=0.5)
        }
        par(xpd=FALSE)
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
