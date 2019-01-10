# Function to sort fits by magnitude of focal parameter point estimates
order.of.fits<-function(ffr.fits, model, parm, order=FALSE){
  if(order){
    foc.parms <- unlist(lapply(ffr.fits, function(x) coef(x$fits[[model]])[parm]))
    how.to.order <- order(foc.parms)
  }else{
    how.to.order <- 1:length(ffr.fits)
  }
}

# Function to plot focal parameter estimates across all fitted datasets
plot.coefs <- function(
                ffr.fits, 
                model=c('Holling.Type.I','Holling.Type.II','Beddington.DeAngelis','Crowley.Martin',
                        'Stouffer.Novak.I','Ratio','Arditi.Ginzburg','Hassell.Varley','Arditi.Akcakaya',
                        'Arditi.Akcakaya.Method.2'),
                parameter=c('attack','handling','interference','phi_denom','exponent'),  # add others later as needed
                ilink=identity,
                plot.SEs=FALSE,
                display.outlier.ests=FALSE,
                xlim=NULL, ... ){

  model <- match.arg(model)
  parameter <- match.arg(parameter)
  ilink <- match.fun(ilink)

  # scaffold of a plot that doesn't actually show anything
  plot(
    y = 1:length(ffr.fits),
    # do we really need this if type='n' anyway?
    # x = unlist(lapply(ffr.fits, function(x){ilink(coef(x$fits[[model]])[parameter])})), 
    x = 1:length(ffr.fits),
    type='n',
    yaxt='n',
    xlim=xlim,
    ...
  )
  
  # tick marks to indicate different data sets
  axis(side=2, at=1:length(ffr.fits), labels=FALSE)
  
  # mark where the existing models fall
  abline(v=0, lty=2)
  abline(v=1, lty=2)
  
  # length of arrows to indicate values off the plot
  delta.arrow <- 0.1
  
  # plot these bad boys
  i <- 1
  for(nn in names(ffr.fits)){
    x <- ffr.fits[[nn]]
    
    # color points depending on predator/parasitoid
    col <- ifelse(x$study.info$predator, "black", "red")	
    
    # the median estimate is easy to determine regardless of the type of data
    mm <- x$estimates[[model]]["50%",parameter,"estimate"]
    
    # cheeky upper and lower bounds in the absence of SE information
    lb <- mm
    ub <- mm
    
    # make all lines the equivalent for now
    lty <- "solid"
    
    # if the point estimate is out of bounds don't even bother trying to profile its uncertainty
    if(mm < xlim[1] | mm > xlim[2]){
      if(mm > xlim[2]){
        arrows(xlim[2]-delta.arrow, i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        if(display.outlier.ests){
            text(xlim[2]-delta.arrow,i,round(mm,1),pos=2,cex=0.7*par()$cex)
        }
      }else{
        arrows(xlim[1]+delta.arrow, i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        if(display.outlier.ests){
          text(xlim[1]+delta.arrow,i,round(mm,1),pos=4,cex=0.7*par()$cex)
        }
      }
    }else{
    
    # Three ways to estimate intervals
    if(plot.SEs){
      
      # if we did not bootstrap then try profile or approximate	
      if(x$estimates[[model]]["n",1,1] == 1){
        
        # (1) estimate the profile confidence interval
            # do so for all model parameters because doing so for focal parameter can cause errors.
        cf <- try(confint(x$fits[[model]],try_harder=TRUE, level=0.68, tol.newmin=Inf, quietly=TRUE))
        
        # if profiling code was successful
        if(!inherits(cf, "try-error")){
          
          # best case is solid line
          lty <- "solid"
          
          lb <- cf[parameter,1]
          ub <- cf[parameter,2]
          
          # sometimes we profile things but still get NA intervals, so stretch interval(s) to extremes of plot
          lb <- ifelse(is.na(lb), xlim[1], lb)
          ub <- ifelse(is.na(ub), xlim[2], ub)
        }else{
          # (2) if profiling was not successful then assume quadratic approximation
         
          # quadratic approximation is dashed line
          lty <- "dashed"
          
          # get the SEs directly from the model output
          se <- coef(summary(x$fits[[model]]))[parameter,"Std. Error"]
          lb <- ifelse(is.na(se), xlim[1], mm - se)
          ub <- ifelse(is.na(se), xlim[2], mm + se)
        }
      }else{
        # (3) if we bootstrapped then use the quantiles
        
        # bootstrapped is dotted line
        lty <- "dotted"
        
        # use the central interval equivalent to one SD as the bounds
        lb <- x$estimates[[model]]["16%",parameter,"estimate"]
        ub <- x$estimates[[model]]["84%",parameter,"estimate"]
      }
      
      # For any of the above, don't plot off the figure
      lb <- ifelse(lb < xlim[1], xlim[1], lb)
      ub <- ifelse(ub > xlim[2], xlim[2], ub)
    }
    
    if(mm > xlim[1] & mm < xlim[2]){
      # draw the error bars
      segments(ilink(lb), i, ilink(ub), i, col=col, lty=lty)
      
      # arrow up the limiting cases
      if(lb==xlim[1]){
        arrows(xlim[1], i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
      }
      if(ub==xlim[2]){
        arrows(xlim[2], i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
      }
      
      # plot the actual mean estimate
      points(y=c(i),x=c(ilink(mm)),col=col,bg=col,pch=21)
    }
    }
    i <- i + 1
  }
}