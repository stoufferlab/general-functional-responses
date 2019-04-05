# Function to sort fits by magnitude of focal parameter point estimates
order.of.fits<-function(ffr.fits, order=FALSE, model=NULL, order.parm=NULL, point.est=c('median','mean')){
  point.est <- match.arg(point.est)
  point.est <- ifelse(point.est=='median', '50%', point.est)
  if(order){
    if(order.parm=='Sample size'){
      foc.parms <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
      how.to.order <- order(foc.parms)
    }else{
        foc.parms <- unlist(lapply(ffr.fits, function(x) x$estimates[[model]][point.est, order.parm, "estimate"]))
        how.to.order <- order(foc.parms)
      }
    }else{
      how.to.order <- 1:length(ffr.fits)
    }
}

# Function to plot focal parameter estimates across all fitted datasets
plot.coefs <- function(
                ffr.fits, 
                model=c('Holling.Type.I',
                        'Holling.Type.II',
                        'Beddington.DeAngelis',
                        'Crowley.Martin',
                        'Stouffer.Novak.I',
                        'Ratio',
                        'Arditi.Ginzburg',
                        'Hassell.Varley',
                        'Arditi.Akcakaya',
                        'Arditi.Akcakaya.Method.2'),
                parameter=c('attack',
                            'handling',
                            'interference',
                            'phi_denom',
                            'exponent'),  # add others later as needed
                ilink=identity,
                point.est=c('median','mean'),
                color.factor=c('None','Parasitoids','Replacement'),
                color.vector=NULL,
                pch.factor=c('None','Parasitoids','Replacement'),
                pch.vector=NULL,
                plot.SEs=FALSE,
                display.outlier.ests=FALSE,
                xlim=NULL,
                labels=FALSE,
                vertLines=c(0,1),
                ... ){

  model <- match.arg(model)
  parameter <- match.arg(parameter)
  ilink <- match.fun(ilink)
  point.est <- match.arg(point.est)
    point.est <- ifelse(point.est=='median', '50%', point.est)
  color.factor <- match.arg(color.factor)
  color.vector <- color.vector
  pch.factor <- match.arg(pch.factor)
  pch.vector <- pch.vector

  if(is.null(xlim)){
    ests <- unlist(lapply(ffr.fits, function(x) x$estimates[[model]][point.est, parameter, "estimate"]))
    xlim <- range(ilink(ests))
  }
  
  # scaffold of a plot that doesn't actually show anything
  plot(
    y = 1:length(ffr.fits),
    x = 1:length(ffr.fits),
    type='n',
    yaxt='n',
    xlim=xlim,
    ...
  )
  
  # tick marks to indicate different data sets
  sample.sizes <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))
  if(labels){
    labels <- names(sample.sizes)
    labels<-sub('./Dataset_Code/','',labels)
    labels<-sub('.{2}$','',labels)
    labels <- paste0(labels, ' (',sample.sizes,')')
  }
  axis(side=2, at=1:length(ffr.fits), labels=labels, cex.axis=0.5)
  
  # mark where the existing models fall
  abline(v=vertLines, lty=2, col='grey')

  # length of arrows to indicate values off the plot
  delta.arrow <- 0.1
  
  # create a progress bar that shows how far along the bootstrapping is
  # require(progress)
  # pb <- progress_bar$new(
  #   format = "  bootstrapping [:bar] :percent eta: :eta",
  #   total = boot.reps,
  #   show_after = 0,
  #   force = TRUE,
  #   clear = FALSE
  # )
  
  # plot these bad boys 
  # (Note: This determines intervals on estimated scale before plotting on backtransformed scale)
  i <- 1
  for(i in 1:length(ffr.fits)){
    nn <- names(ffr.fits)[i]
    x <- ffr.fits[[nn]]
    
    # color points depending on factor or vector of choice
    col <- 'black'
    if(!color.factor=='None'){
      if(color.factor=='Parasitoids' & !x$study.info$predator){
        col <- 'red'  }
      if(color.factor=='Replacement' & x$study.info$replacement){
        col <- 'red'  }
    }
    if(length(color.vector)==length(ffr.fits)){  # note that this overrides color.factor
      col <- color.vector[i]
    }
    
    # set pch of points depending on factor or vector of choice
    pch <- 19
    bg <- col
    if(!pch.factor=='None'){
      if(pch.factor=='Parasitoids' & !x$study.info$predator){
        pch <- 21
        bg <- 'white'
        }
      if(pch.factor=='Replacement' & x$study.info$replacement){
        pch <- 21
        bg <- 'white'
        }
    }
    if(length(pch.vector)==length(ffr.fits)){  # note that this overrides pchor.factor
      pch <- pch.vector[i]
    }

    # make all lines the equivalent for now
    lty <- "solid"
    
    # extract point estimate (the median estimate is easy to determine regardless of the type of data so should be default)
    mm <- x$estimates[[model]][point.est, parameter, "estimate"]
    mm.link <- ilink(mm)
    
    # cheeky upper and lower bounds in the absence of SE information
    lb <- mm
    ub <- mm

    # if the point estimate is out of bounds don't even bother trying to profile its uncertainty
    if(mm.link < xlim[1] | mm.link > xlim[2]){
      if(mm.link > xlim[2]){
        arrows(xlim[2]-delta.arrow, i, xlim[2]+delta.arrow, i,length=delta.arrow*0.66, col=col, lty=lty)
        if(display.outlier.ests){
            text(xlim[2]-delta.arrow,i,round(mm.link,1), pos=2,cex=0.7*par()$cex)
        }
      }else{
        arrows(xlim[1]+delta.arrow, i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        if(display.outlier.ests){
            text(xlim[1]+delta.arrow,i,round(mm.link,1), pos=4,cex=0.7*par()$cex)
        }
      }
    }else{
    
    # Three ways to estimate intervals
    if(plot.SEs){
        # if we did not bootstrap then try (1) profile or (2) approximate	
        if(x$estimates[[model]]["n",1,1] == 1){ 
          
          # (1) estimate the profile confidence interval
              # do so for all model parameters because doing so for focal parameter can cause errors
          if(model!='Arditi.Akcakaya.Method.2'){
             cf <- try(confint(x$fits[[model]], try_harder=TRUE, level=0.68, tol.newmin=Inf, quietly=TRUE))
          }else{
            cf <- TRUE
            class(cf) <- 'try-error'
          }
          
          # if profiling code was successful
          if(!inherits(cf, "try-error")){
            
            # best case is solid line
            lty <- "solid"
            
            lb <- cf[parameter,1]
            ub <- cf[parameter,2]
  
          }else{
          # (2) if profiling is unsuccessful or AA2 method was used then assume quadratic approximation

            # quadratic approximation is dashed line
            lty <- "dashed"
            
            # get the SEs directly from the model output
            if(model=='Arditi.Akcakaya.Method.2'){
              se <- x$estimates$Arditi.Akcakaya.Method.2[,,'std.error'][point.est,parameter]
            }else{
              se <- coef(summary(x$fits[[model]]))[parameter,"Std. Error"]
            }
            lb <- mm - se
            ub <- mm + se
          }
        }else{# (3) if we bootstrapped then use the quantiles
          
          # bootstrapped is dotted line
          lty <- "dotted"
          
          # use the central interval equivalent to one SD as the bounds
          lb <- x$estimates[[model]]["16%",parameter,"estimate"]
          ub <- x$estimates[[model]]["84%",parameter,"estimate"]
        }
        
        lb.link <- ilink(lb)
        ub.link <- ilink(ub)
        
        # sometimes se is NA or we profile things but still get NA intervals, so stretch interval(s) to extremes of plot
        lb.link <- ifelse(is.na(lb.link), xlim[1], lb.link)
        ub.link <- ifelse(is.na(ub.link), xlim[2], ub.link)
        
        # For any of the above, don't plot off the figure
        lb.link <- ifelse(lb.link < xlim[1], xlim[1], lb.link)
        ub.link <- ifelse(ub.link > xlim[2], xlim[2], ub.link)
      
      
      if(mm.link > xlim[1] & mm.link < xlim[2]){
        # draw the error bars
        segments(lb.link, i, ub.link, i, col=col, lty=lty)
        
        # arrow up the limiting cases
        if(lb.link <= xlim[1] & parameter!='exponent'){
          arrows(xlim[1], i, xlim[1]-delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        }
        if(ub.link >= xlim[2]){
          arrows(xlim[2], i, xlim[2]+delta.arrow, i, length=delta.arrow*0.66, col=col, lty=lty)
        }
      }
    }
      # plot the actual estimate
      points(y=i, x=mm.link,
             col=col, pch=pch, bg=bg)
    
    }
    
    # update the progress bar
    # pb$tick()
    
    i <- i + 1
  }
}