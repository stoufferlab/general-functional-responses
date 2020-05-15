
# essential libraries
library(RColorBrewer)
library(shape)

# Function for number reformatting
format.number <-function(number){
  if(number < -100){
    result <- format(signif(number,3), digits=3, nsmall=2, scientific=TRUE)
  }else if(number < -1){
    result <- format(round(number, 2), nsmall=2)
  }else if(number < -0.001){
    result <- format(signif(number,3), nsmall=0, scientific=TRUE)
  }else if(number < 0.001){
    result <- format(signif(number,3), nsmall=0, scientific=TRUE)
  }else if(number < 1){
    result <- format(round(number,3),nsmall = 3)
  }else if(number < 100){
    result <-  format(round(number, 2), nsmall=2)
  }else{ # > 100
    result <- format(signif(number,3) , nsmall=2, scientific=TRUE)
  }
  return(result)
}

# read in the profiled CIs
load(
	file='../../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata'
)

# reverse order
ffr.cfs <- ffr.cfs[rev(1:length(ffr.cfs))]

pdf(
    '../../../../results/R/OnePredTwoPrey_figs/OnePredTwoPrey_phi.pdf',
    height=2.875,
    width=2.25
)

par(
    mar=c(3,7,2.5,0.75),
    oma = c(0, 0, 0, 0),
    mgp=c(1.25,0.1,0.0),
    tcl=-0.1,
    las=1,
    cex=0.7,
    yaxs='i'
)

# x axis limits
xlim <- c(-5,5)

# some colors
CR<-brewer.pal(n = 9, name = 'Blues')

# jitter
dy <- 0.4

# dataset labels
labels <- unlist(lapply(ffr.cfs, function(x) x$study.info$datasetName))
labels<-gsub('_',' ',labels)
labels <- paste(labels,"  ")


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

	# tick marks to indicate different data sets
	axis(side=2, at=1:length(ffr.cfs), labels=labels, cex.axis=0.5, las=1, lwd=0, lwd.ticks=1)
	axis(side=1, cex.axis=0.7, mgp=c(1.25,0,0), lwd=0, lwd.ticks=1)
	box(lwd=1)

	# length of arrows to indicate values off the plot
	delta.arrow <- 0.04*diff(xlim)

	# plot these bad boys
	for(i in 1:length(ffr.cfs)){
		x <- ffr.cfs[[i]]
		col <- ifelse(
			any(rank(x$AICs)[c(6,7)]==1),
			"black",
			CR[7]
		)
		for(j in 1:nrow(x$profile)){
			# the median estimate is easy to determine regardless of the type of data
			mm <- x$profile[j,"est"]
			lb <- x$profile[j,"lb"]
			ub <- x$profile[j,"ub"]

			# make all lines the equivalent for now
			lty <- switch(x$profile[j,"method"],
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
					cex=0.33
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
					cex=0.33
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
				# 	arrows(xlim[2]-delta.arrow, i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
				# }else{
				# 	arrows(xlim[1]+delta.arrow, i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
				# }
			# }
		}
	}
# }

# xlim <- c(-5, 5)

# modeltype <- "Holling II Hybrid Hybrid"
# parameters <- c("phi_ij","phi_ji")

# # dev.new()

# plot.coefs(
# 	ffr.fits,
# 	modeltype = modeltype,
# 	parameters = parameters,
# 	plot.SEs = TRUE,
# 	xlim = xlim,
# 	ilink=identity,
# 	xlab="Effect of feeding on feeding",
# 	ylab="Dataset"
# )

# dev.copy2pdf(file="../../../results/R/OnePredTwoPrey_Phis.pdf")
dev.off()
