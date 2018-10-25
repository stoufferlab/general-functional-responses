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
source('../lib/AA_method.R')
source('../lib/holling_method_one_predator_one_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEBUG: we should probably move these functions elsewhere so this file only calls functions/runs analysis
# library(devtools)
# devtools::install_github("bbolker/broom")
library(broom) # for tidy()

mytidy <- function(fit){ 
  tfit <- tidy(fit)
  terms <- tfit$term
  tfit$term <- NULL
  out <- matrix(as.numeric(unlist(tfit)), nrow=nrow(tfit), ncol=ncol(tfit), byrow=FALSE)
  rownames(out) <- terms
  colnames(out) <- colnames(tfit)
  out
}

make.array <- function(ffr.fit,boot.reps){
  t.ffr.fit <- mytidy(ffr.fit)
  out <- array(NA, dim=c(dim(t.ffr.fit), boot.reps))
  dimnames(out) <- dimnames(t.ffr.fit)
  out
}

summarize.boots <- function(x){
  c(mean=mean(x, na.rm=TRUE), quantile(x,c(0.025,0.16,0.5,0.84,0.975),na.rm=TRUE), n=sum(!is.na(x)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# master list of datasets
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)

# remove template files which don't actually read data
datasets <- grep("template",datasets,invert=TRUE,value=TRUE)

# remove zzz files which are placeholders while a dataset is being cleaned/incorporated
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)

# create a container for the things that get fit
ffr.fits <- list()

# # DEBUG: for testing only
# datasets <- c("./Dataset_Code/Prokopenko_2017.R")

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

	if(!grepl("H", this.study$runswith)){
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

		library(progress)
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

	  	  	# DEBUG sometimes the fits fail for bootstrapped data; we could just skip that replicate and continue on?
	    	# fit a series of functional response models
	    	success <- try({
	    		ffr.hollingI <- fit.holling.like(d, this.study, "Holling I")
				ffr.hollingII <- fit.holling.like(d, this.study, "Holling II")
				ffr.bd <- fit.holling.like(d, this.study, "Beddington-DeAngelis")
				ffr.cm <- fit.holling.like(d, this.study, "Crowley-Martin")
				ffr.sn1 <- fit.holling.like(d, this.study, "Stouffer-Novak I")
	    		# ifelse(okay4AAmethod(d), fit.AAmethod <- AAmethod(d,expttype), fit.AAmethod <- NA)
	    	})
	    	
	    	if(!inherits(success, "try-error")){
		    	if(b == 1){
		    		# Containers for parameter estimates
					boots.HT.I <- make.array(ffr.hollingI, boot.reps)
					boots.HT.II <- make.array(ffr.hollingII, boot.reps)
					boots.BD <- make.array(ffr.bd, boot.reps)
					boots.CM <- make.array(ffr.cm, boot.reps)
					boots.SN.I <- make.array(ffr.sn1, boot.reps)
		    	}

		    	# add the parameters for this fit to the array of fits
				boots.HT.I[,,b] <- mytidy(ffr.hollingI)
				boots.HT.II[,,b] <- mytidy(ffr.hollingII)
				boots.BD[,,b] <- mytidy(ffr.bd)
				boots.CM[,,b] <- mytidy(ffr.cm)
				boots.SN.I[,,b] <- mytidy(ffr.sn1)
				# boots.SN.Numer[,,b] <- mytidy(ffr.sn2)
				# boots.SN.III[,,b] <- mytidy(ffr.sn3)
				# boots.HV[,,b] <- mytidy(ffr.hv)
				# boots.AG[,,b] <- mytidy(ffr.ag)
				# boots.AA[,,b] <- mytidy(ffr.aa)
				# if(okay4AAmethod(d)){ boots.AA2[,,b] <- fit.AAmethod$estimates  }

				pb$tick()
				b <- b + 1
			}
	  	}
	  
		# ~~~~~~~~~~~~~~~~~~~~
		# Summarize bootstraps
		# ~~~~~~~~~~~~~~~~~~~~
		HT.I.ests <- as.array(apply(boots.HT.I, c(1,2), summarize.boots))
		HT.II.ests <- as.array(apply(boots.HT.II, c(1,2), summarize.boots))
		BD.ests <- as.array(apply(boots.BD, c(1,2), summarize.boots))
		CM.ests <- as.array(apply(boots.CM, c(1,2), summarize.boots))
		SN.I.ests <- as.array(apply(boots.SN.I, c(1,2), summarize.boots))
		# SN.Numer.ests <- as.array(apply(boots.SN.Numer,c(1,2), summarize.boots))
		# # SN.III.ests <- as.array(apply(boots.SN.III,c(1,2), summarize.boots))
		# HV.ests <- as.array(apply(boots.HV,c(1,2), summarize.boots))	  
		# AG.ests <- as.array(apply(boots.AG,c(1,2), summarize.boots))
		# AA.ests <- as.array(apply(boots.AA,c(1,2), summarize.boots))
		# if(okay4AAmethod(d)){ 
		#     AA2.ests <- as.array(apply(boots.AA2,c(1,2), summarize.boots))
		# }else{
		# 	AA2.ests <- NA
		# }
  
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
				Stouffer.Novak.I = ffr.sn1
	          	# Stouffer.Novak.Numer = ffr.sn2,
	          	# Stouffer.Novak.III = ffr.sn3,
	          	# Hassell.Varley = ffr.hv,
	          	# Arditi.Ginzburg = ffr.ag,
	          	# Arditi.Akcakaya = ffr.aa,
	          	# Arditi.Akcakaya.Method.2 = fit.AAmethod
	        ),
	        # boots = list(
			# 	Holling.TypeI = boots.HT.I,
			# 	Holling.TypeII = boots.HT.II,
			# 	Beddington.DeAngelis = boots.BD,
			# 	Crowley.Martin = boots.CM,
			# 	Stouffer.Novak.I = boots.SN.I
			# ),
			estimates = list(
			    Holling.Type.I = HT.I.ests,
			    Holling.Type.II = HT.II.ests,
			    Beddington.DeAngelis = BD.ests,
			    Crowley.Martin = CM.ests,
			    Stouffer.Novak.I = SN.I.ests
			    # Stouffer.Novak.Numer = SN.Numer.ests,
			    # Stouffer.Novak.III = SN.III.ests,
			    # Hassell.Varley = HV.ests,
			    # AArditi.Ginzburg = G.ests,
			    # Arditi.Akcakaya = AA.ests,
			    # Arditi.Akcakaya.Method.2 = AA2.ests)
			)
		)
	# })
	}

	# break
}

# save the mega container which includes all FR fits
save(ffr.fits,file='../../../results/R/ffr.fits_OnePredOnePrey.Rdata')
