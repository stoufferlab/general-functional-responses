rm(list = ls())
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
datasets <- c("./Dataset_Code/Mertz_1968.R")


# create a container for the things that get fit
ffr.fits <- list()

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
	# NOTE: optimization is on log-transformed values
	#############################################	 

	# H: holling-like, R: ratio-like, T: test set (or combinations thereof)
	if(!grepl("RT", this.study$runswith)){ 
	# if(!grepl("H|R", this.study$runswith)){ 
		message(paste0("No to ",datasets[i]))
	}else{
		# print out which dataset is being analyzed
		message(paste0("Yes to ",datasets[i]))

		# Do data need to be bootstrapped?
		if("Nconsumed.mean" %in% colnames(d)){
			boot.reps <- 3# 250
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
		    		# creates containers for parameter estimates using first rep as a template
		    		boots.HT.I <- boots.HT.II <- boots.BD <- boots.CM <- boots.SN.I <- array(NA, c(1,1,boot.reps))
		    		if(grepl("H", this.study$runswith)){
	  					boots.HT.I <- make.array(ffr.hollingI, boot.reps)
	  					boots.HT.II <- make.array(ffr.hollingII, boot.reps)
	  					boots.BD <- make.array(ffr.bd, boot.reps)
	  					boots.CM <- make.array(ffr.cm, boot.reps)
	  					boots.SN.I <- make.array(ffr.sn1, boot.reps)
	  				}
  					
  					boots.R <- boots.AG <- boots.HV <- boots.AA <- boots.AAm <- array(NA,c(1,1,boot.reps))
    				if(grepl("R", this.study$runswith)){
      					boots.R <- make.array(ffr.ratio, boot.reps)
      					boots.AG <- make.array(ffr.ag, boot.reps)
      					boots.HV <- make.array(ffr.hv, boot.reps)
      					boots.AA <- make.array(ffr.aa, boot.reps)
      					if(okay4AAmethod(d)){
      					  boots.AAm <- make.array(ffr.aam$estimates, boot.reps)
    					  }
    				}
		      	}

		    	# add the parameters for this rep to the array of fits
		    	boots.HT.I[,,b] <- boots.HT.II[,,b] <- boots.BD[,,b] <- boots.CM[,,b] <- boots.SN.I[,,b] <- NA #array(NA, c(1,1))
				if(grepl("H", this.study$runswith)){
					boots.HT.I[,,b] <- mytidy(ffr.hollingI)
					boots.HT.II[,,b] <- mytidy(ffr.hollingII)
					boots.BD[,,b] <- mytidy(ffr.bd)
					boots.CM[,,b] <- mytidy(ffr.cm)
					boots.SN.I[,,b] <- mytidy(ffr.sn1)
					# boots.SN.Numer[,,b] <- mytidy(ffr.sn2)
					# boots.SN.III[,,b] <- mytidy(ffr.sn3)
				}
				
				boots.R[,,b] <- boots.AG[,,b] <- boots.HV[,,b] <- boots.AA[,,b] <- boots.AAm[,,b] <- NA #array(NA, c(1,1))
				if(grepl("R", this.study$runswith)){
	  				boots.R[,,b] <- mytidy(ffr.ratio)
	  				boots.AG[,,b] <- mytidy(ffr.ag)
	  				boots.HV[,,b] <- mytidy(ffr.hv)
	  				boots.AA[,,b] <- mytidy(ffr.aa)
	  				if(okay4AAmethod(d)){ 
	  				    boots.AAm[,,b] <- ffr.aam$estimate
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
		HT.I.ests <- HT.II.ests <- BD.ests <- CM.ests <- SN.I.ests <- NA
		if(grepl("H", this.study$runswith)){
			HT.I.ests <- as.array(apply(boots.HT.I, c(1,2), summarize.boots))
			HT.II.ests <- as.array(apply(boots.HT.II, c(1,2), summarize.boots))
			BD.ests <- as.array(apply(boots.BD, c(1,2), summarize.boots))
			CM.ests <- as.array(apply(boots.CM, c(1,2), summarize.boots))
			SN.I.ests <- as.array(apply(boots.SN.I, c(1,2), summarize.boots))
			# SN.Numer.ests <- as.array(apply(boots.SN.Numer,c(1,2), summarize.boots))
			# SN.III.ests <- as.array(apply(boots.SN.III,c(1,2), summarize.boots))
		}
		
		R.ests <- AG.ests <- HV.ests <- AA.ests <- AAm.ests <- NA
		if(grepl("R", this.study$runswith)){
	  		R.ests <- as.array(apply(boots.R,c(1,2), summarize.boots))
	  		AG.ests <- as.array(apply(boots.AG,c(1,2), summarize.boots))
	  		HV.ests <- as.array(apply(boots.HV,c(1,2), summarize.boots))
	  		AA.ests <- as.array(apply(boots.AA,c(1,2), summarize.boots))
        if(okay4AAmethod(d)){
            AAm.ests <- as.array(apply(boots.AAm,c(1,2), summarize.boots))
        }
		}
	  	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# save the (last) fits and some data aspects
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ffr.fits[[datasets[i]]] <- list(
	  		study.info = c(
	  			datadir = datadir,
	  			sample.size = nrow(d),
	  			data=d,
	  	        this.study  	        
	  	    ),
			fits = c(
  				Holling.Type.I = ffr.hollingI,
  				Holling.Type.II = ffr.hollingII,
  				Beddington.DeAngelis = ffr.bd,
  				Crowley.Martin = ffr.cm,
  				Stouffer.Novak.I = ffr.sn1,
        		# Stouffer.Novak.Numer = ffr.sn2,
        		# Stouffer.Novak.III = ffr.sn3,
				Ratio = ffr.ratio,
				Arditi.Ginzburg = ffr.ag,
				Hassell.Varley = ffr.hv,
				Arditi.Akcakaya = ffr.aa,
				Arditi.Akcakaya.Method.2 = ffr.aam
			),
      boots = list(
      	Holling.TypeI = boots.HT.I,
      	Holling.TypeII = boots.HT.II,
      	Beddington.DeAngelis = boots.BD,
      	Crowley.Martin = boots.CM,
      	Stouffer.Novak.I = boots.SN.I,
      	Ratio = boots.R,
        Arditi.Ginzburg = boots.AG,
        Hassell.Varley = boots.HV,
        Arditi.Akcakaya = boots.AA,
        Arditi.Akcakaya.Method.2 = boots.AAm
        ),
			estimates = list(
			    Holling.Type.I = HT.I.ests,
			    Holling.Type.II = HT.II.ests,
			    Beddington.DeAngelis = BD.ests,
			    Crowley.Martin = CM.ests,
			    Stouffer.Novak.I = SN.I.ests,
			    # Stouffer.Novak.Numer = SN.Numer.ests,
			    # Stouffer.Novak.III = SN.III.ests,
			    Ratio = R.ests,
			    Arditi.Ginzburg = AG.ests,
			    Hassell.Varley = HV.ests,
			    Arditi.Akcakaya = AA.ests,
			    Arditi.Akcakaya.Method.2 = AAm.ests
			)
		)
	# })
	}

	# break
}

# save the mega container which includes all FR fits
save(ffr.fits,file='../../../results/R/ffr.fits_OnePredOnePrey.Rdata')
