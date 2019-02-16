rm(list = ls())
sinkMessages <- TRUE # set to FALSE if you want to match messages in real time 
# or TRUE to have them silently saved to file instead.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
	Sys.getenv("LOGNAME"),
	stouffer = '../../../dropbox_data/Data',
	marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a few utility functions
source('../lib/study_info.R')
source('../lib/bootstrap_data.R')
source('../lib/mytidySumm.R')
source('../lib/AA_method.R')
source('../lib/holling_method_one_predator_one_prey.R')
source('../lib/ratio_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(progress)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# master list of datasets
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)

# remove template files which don't actually read data
datasets <- grep("template",datasets,invert=TRUE,value=TRUE)

# remove zzz files which are placeholders while a dataset is being cleaned/incorporated
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)

# # DEBUG: for testing only
datasets <- c("./Dataset_Code/Walde_1984.R")  # Occasional Hessian problem


# create mega container for the things that get fit 
ffr.ests <- ffr.fits <- list()

# start capturing the progress and warning messages
if(sinkMessages){
  options(warn=1)
  Mesgs <- file('../../../results/R/ffr.fits_OnePredOnePrey_LOG.txt', open='wt')
  sink(Mesgs, type="message")
}

# fit everything on a dataset by dataset basis
for(i in 1:length(datasets)){
	# loads the data into data frame 'd' and specifies data-specific parameters
	source(datasets[i])

	# grab some info from the google doc
	this.study <- study.info(datadir)

	# put all datasets into terms of hours
	if(!is.null(d$Time)){
		d$Time <- switch(this.study$timeunits,
			Seconds = d$Time / 3600.,
			Minutes = d$Time / 60.,
			Hours = d$Time,
			Days = d$Time * 24,
			Unavailable = rep(NA, nrow(d))
		)
	}

	# save original data in case we need to bootstrap it
	d.orig <- d

	#############################################
	# fit all the functional response models
	# NOTE: optimization is on log-transformed values for some of the parameters
	#############################################	 

	# H: holling-like, R: ratio-like, T: test set (or combinations thereof)
	# if(!grepl("R", this.study$runswith)){ 
	if(!grepl("H|R", this.study$runswith)){
		message(paste0("No to ",datasets[i]))
	}else{
		# print out which dataset is being analyzed
		message(paste0("Yes to ",datasets[i]))

		# Do data need to be bootstrapped?
		if("Nconsumed.mean" %in% colnames(d)){
			boot.reps <- 250
		}else{
			boot.reps <- 1
		}

		# create a progress bar that shows how far along the bootstrapping is
		pb <- progress_bar$new(
			format = "  bootstrapping [:bar] :percent eta: :eta",
			total = boot.reps,
			show_after = 0,
			force = TRUE,
			clear = FALSE
		)

		# perform boot.reps fits depending on the nature of the data
		b <- 1
	  	while(b <= boot.reps){
	  		# generate bootstrapped data if necessary
			if("Nconsumed.mean" %in% colnames(d.orig)){
	  			d <- bootstrap.data(d.orig, this.study$replacement)
	  	  	}

	  	  # fit a series of functional response models
	  	  # (sometimes the fits fail for bootstrapped data. we skip that replicate and continue on)
	    	success <- try({
	    		ffr.hollingI <- ffr.hollingII <- ffr.bd <- ffr.cm <- ffr.sn1 <- array(NA,c(1,1))
	    		if(grepl("H", this.study$runswith)){
	  	    	ffr.hollingI <- fit.holling.like(d, this.study, "Holling I")
	  				ffr.hollingII <- fit.holling.like(d, this.study, "Holling II")
	  				ffr.bd <- fit.holling.like(d, this.study, "Beddington-DeAngelis")
	  				ffr.cm <- fit.holling.like(d, this.study, "Crowley-Martin")
	  				ffr.sn1 <- fit.holling.like(d, this.study, "Stouffer-Novak I")
	  			}
  				
  				ffr.ratio <- ffr.ag <- ffr.hv <- ffr.aa <- ffr.aam <- array(NA,c(1,1))
  				if(grepl("R", this.study$runswith)){
    				ffr.ratio <- fit.ratio.like(d, this.study, "Ratio")
    				ffr.ag <- fit.ratio.like(d, this.study, "Arditi-Ginzburg")
    				ffr.hv <- fit.ratio.like(d, this.study, "Hassell-Varley")
    				ffr.aa <- fit.ratio.like(d, this.study, "Arditi-Akcakaya")
  	    		if(okay4AAmethod(d)){
  	    		   ffr.aam <- AAmethod(d,this.study$replacement)
  	    		}
  				}
	    	})
	    	
	    	if(!inherits(success, "try-error")){
		    	if(b == 1){
		    	  # create containers for AIC values of all bootstrap fits
		    	  aic.hollingI <- aic.hollingII <- aic.bd <- aic.cm <- aic.sn1 <- vector("numeric",boot.reps)
		    	  # aic.sn2 <- aic.sn3 <- vector("list",boot.reps)
		    	  aic.ratio <- aic.ag <- aic.hv <- aic.aa <- vector("numeric",boot.reps)
		    	  
		    	  # create containers for parameter estimates using first rep as a template
		    		boots.hollingI <- boots.hollingII <- boots.bd <- array(NA, c(1,1,boot.reps))
		    		boots.cm <- boots.sn1 <- array(NA, c(1,1,boot.reps))
		    		# boots.sn2 <- boots.sn3 <- array(NA, c(1,1,boot.reps))
		    		
		    		if(grepl("H", this.study$runswith)){
	  					boots.hollingI <- make.array(ffr.hollingI, boot.reps)
	  					boots.hollingII <- make.array(ffr.hollingII, boot.reps)
	  					boots.bd <- make.array(ffr.bd, boot.reps)
	  					boots.cm <- make.array(ffr.cm, boot.reps)
	  					boots.sn1 <- make.array(ffr.sn1, boot.reps)
	  					# boots.sn2 <- make.array(ffr.sn2, boot.reps)
	  					# boots.sn3 <- make.array(ffr.sn3, boot.reps)
	  				}
  					
  					boots.ratio <- boots.ag <- boots.hv <- boots.aa <- boots.aam <- array(NA,c(1,1,boot.reps))
    				if(grepl("R", this.study$runswith)){
      					boots.ratio <- make.array(ffr.ratio, boot.reps)
      					boots.ag <- make.array(ffr.ag, boot.reps)
      					boots.hv <- make.array(ffr.hv, boot.reps)
      					boots.aa <- make.array(ffr.aa, boot.reps)
      					if(okay4AAmethod(d)){
      					  boots.aam <- make.array(ffr.aam$estimates, boot.reps)
    					  }
    				}
		      	}

	    	  # Save the AIC values for the rep
	    	  if(grepl("H", this.study$runswith)){
	    	    aic.hollingI[[b]] <- AIC(ffr.hollingI)
	    	    aic.hollingII[[b]] <- AIC(ffr.hollingII)
	    	    aic.bd[[b]] <- AIC(ffr.bd)
	    	    aic.cm[[b]] <- AIC(ffr.cm)
	    	    aic.sn1[[b]] <- AIC(ffr.sn1)
	    	    # aic.sn2[[b]] <- AIC(ffr.sn2)
	    	    # aic.sn3[[b]] <- AIC(ffr.sn3)
	    	  }
	    	  if(grepl("R", this.study$runswith)){
	    	    aic.ratio[[b]] <- AIC(ffr.ratio)
	    	    aic.ag[[b]] <- AIC(ffr.ag)
	    	    aic.hv[[b]] <- AIC(ffr.hv)
	    	    aic.aa[[b]] <- AIC(ffr.aa)
	    	  }
	    	  
	    	  
		    	# Save the estimates for this rep to the array of estimates
		    	boots.hollingI[,,b] <- boots.hollingII[,,b] <- boots.bd[,,b] <- NA
		    	boots.cm[,,b] <- boots.sn1[,,b] <- NA
		    	# boots.sn2[,,b] <- boots.sn2[,,b] <- NA
		    	
  				if(grepl("H", this.study$runswith)){
  					boots.hollingI[,,b] <- mytidy(ffr.hollingI)
  					boots.hollingII[,,b] <- mytidy(ffr.hollingII)
  					boots.bd[,,b] <- mytidy(ffr.bd)
  					boots.cm[,,b] <- mytidy(ffr.cm)
  					boots.sn1[,,b] <- mytidy(ffr.sn1)
  					# boots.SN.Numer[,,b] <- mytidy(ffr.sn2)
  					# boots.SN.III[,,b] <- mytidy(ffr.sn3)
  				}
				
  				boots.ratio[,,b] <- boots.ag[,,b] <- boots.hv[,,b] <- boots.aa[,,b] <- boots.aam[,,b] <- NA 
  				if(grepl("R", this.study$runswith)){
  	  				boots.ratio[,,b] <- mytidy(ffr.ratio)
  	  				boots.ag[,,b] <- mytidy(ffr.ag)
  	  				boots.hv[,,b] <- mytidy(ffr.hv)
  	  				boots.aa[,,b] <- mytidy(ffr.aa)
  	  				if(okay4AAmethod(d)){ 
  	  				    boots.aam[,,b] <- ffr.aam$estimate
  	  				}
  				}

				# update the progress bar
				pb$tick()
				b <- b + 1
			}
	  	}
	  
		# ~~~~~~~~~~~~~~~~~~~~
		# Summarize bootstraps
		# ~~~~~~~~~~~~~~~~~~~~
		ests.hollingI <- ests.hollingII <- ests.bd <- ests.cm <- ests.sn1 <- NA
		if(grepl("H", this.study$runswith)){
			ests.hollingI <- as.array(apply(boots.hollingI, c(1,2), summarize.boots))
			ests.hollingII <- as.array(apply(boots.hollingII, c(1,2), summarize.boots))
			ests.bd <- as.array(apply(boots.bd, c(1,2), summarize.boots))
			ests.cm <- as.array(apply(boots.cm, c(1,2), summarize.boots))
			ests.sn1 <- as.array(apply(boots.sn1, c(1,2), summarize.boots))
			# ests.sn2 <- as.array(apply(boots.sn2,c(1,2), summarize.boots))
			# ests.sn3 <- as.array(apply(boots.sn3,c(1,2), summarize.boots))
			
			AIC.hollingI <- summarize.boots(aic.hollingI)
			AIC.hollingII <- summarize.boots(aic.hollingI)
			AIC.bd <- summarize.boots(aic.bd)
			AIC.cm <- summarize.boots(aic.cm)
			AIC.sn1 <- summarize.boots(aic.sn1)
			# AIC.sn2 <- summarize.boots(aic.sn2)
			# AIC.sn3 <- summarize.boots(aic.sn3)
		}
		
	  ests.ratio <- ests.ag <- ests.hv <- ests.aa <- ests.aam <- NA
		if(grepl("R", this.study$runswith)){
	  	  ests.ratio <- as.array(apply(boots.ratio,c(1,2), summarize.boots))
	  		ests.ag <- as.array(apply(boots.ag,c(1,2), summarize.boots))
	  		ests.hv <- as.array(apply(boots.hv,c(1,2), summarize.boots))
	  		ests.aa <- as.array(apply(boots.aa,c(1,2), summarize.boots))
        if(okay4AAmethod(d)){
            ests.aam <- as.array(apply(boots.aam,c(1,2), summarize.boots))
        }
	  		
	  		AIC.ratio <- summarize.boots(aic.ratio)
	  		AIC.ag <- summarize.boots(aic.ag)
	  		AIC.hv <- summarize.boots(aic.hv)
	  		AIC.aa <- summarize.boots(aic.aa)
		}
		
	  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# save study info, the (last) fits, bootstrap estimates and estimate summaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ffr.fits[[datasets[i]]] <- list(
	  	study.info = c(
            	  			datadir = datadir,
            	  			sample.size = nrow(d.orig),
            	  			data=d.orig,
            	  	    this.study  	        
          	  	    ),
	  	fits = list(
          	  	  Holling.Type.I = ffr.hollingI,
          	  	  Holling.Type.II = ffr.hollingII,
          	  	  Beddington.DeAngelis = ffr.bd,
          	  	  Crowley.Martin = ffr.cm,
          	  	  Stouffer.Novak.I = ffr.sn1,
          	  	  # Stouffer.Novak.II = ffr.sn2,
          	  	  # Stouffer.Novak.III = ffr.sn3,
          	  	  Ratio = ffr.ratio,
          	  	  Arditi.Ginzburg = ffr.ag,
          	  	  Hassell.Varley = ffr.hv,
          	  	  Arditi.Akcakaya = ffr.aa
          	  	),
      boots = list(
                  	Holling.Type.I = boots.hollingI,
                  	Holling.Type.II = boots.hollingII,
                  	Beddington.DeAngelis = boots.bd,
                  	Crowley.Martin = boots.cm,
                  	Stouffer.Novak.I = boots.sn1,
                  	# Stouffer.Novak.II = boots.sn2,
                  	# Stouffer.Novak.III = boots.sn3,
                  	Ratio = boots.ratio,
                    Arditi.Ginzburg = boots.ag,
                    Hassell.Varley = boots.hv,
                    Arditi.Akcakaya = boots.aa,
                    Arditi.Akcakaya.Method.2 = boots.aam
                  ),
	  	AIC = list(
            	  	  Holling.Type.I = AIC.hollingI,
            	  	  Holling.Type.II = AIC.hollingII,
            	  	  Beddington.DeAngelis = AIC.bd,
            	  	  Crowley.Martin = AIC.cm,
            	  	  Stouffer.Novak.I = AIC.sn1,
            	  	  # Stouffer.Novak.II = AIC.sn2,
            	  	  # Stouffer.Novak.III = AIC.sn3,
            	  	  Ratio = AIC.ratio,
            	  	  Arditi.Ginzburg = AIC.ag,
            	  	  Hassell.Varley = AIC.hv,
            	  	  Arditi.Akcakaya = AIC.aa
	  	          ),
			estimates = list(
              			    Holling.Type.I = ests.hollingI,
              			    Holling.Type.II = ests.hollingII,
              			    Beddington.DeAngelis = ests.bd,
              			    Crowley.Martin = ests.cm,
              			    Stouffer.Novak.I = ests.sn1,
              			    # Stouffer.Novak.Numer = ests.sn2,
              			    # Stouffer.Novak.III = ests.sn3,
              			    Ratio = ests.ratio,
              			    Arditi.Ginzburg = ests.ag,
              			    Hassell.Varley = ests.hv,
              			    Arditi.Akcakaya = ests.aa,
              			    Arditi.Akcakaya.Method.2 = ests.aam
            		    	)
	          )
	}

}

# save the mega container
save(ffr.fits,file='../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

if(sinkMessages){
  sink(type="message")
  close(Mesgs)
  options(warn=0)
  readLines("../../../results/R/ffr.fits_OnePredOnePrey_LOG.txt")
}


