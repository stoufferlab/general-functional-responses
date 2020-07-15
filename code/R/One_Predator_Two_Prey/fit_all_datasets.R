rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
	Sys.getenv("LOGNAME"),
	stouffer = '~/Dropbox/Projects/GenFuncResp/Data',
	marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a few utility functions
source('../lib/read_data.R')
source('../lib/study_info.R')
source('../lib/bootstrap_data.R')
source('../lib/mytidySumm.R')
source('../lib/plot_coefs.R')
source('../lib/holling_method_one_predator_two_prey.R')
source('../lib/resid_metrics.R')
# source('../lib/RMSD.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################
library(doParallel)
registerDoParallel(cores=6)
####################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# master list of datasets
datasets <- list.files('./Dataset_Code', pattern=".R$", full.names=TRUE, include.dirs=FALSE)

# define the models which are to be fit
which.models <- c(
	"Holling I",
	"Holling II Specialist Specialist",
	# "Holling II Specialist Hybrid",
	# "Holling II Hybrid Specialist",
	"Holling II Specialist Generalist",
	"Holling II Generalist Specialist",
	"Holling II Generalist Generalist",
	# "Holling II Generalist Hybrid",
	# "Holling II Hybrid Generalist",
	"Holling II Hybrid Hybrid"
)

# fit everything on a dataset by dataset basis
for(i in seq_along(datasets)){
	# create a short nickname for the dataset
	datasetsName <- sub('*.R$','', sub('*./Dataset_Code/','', datasets[i]))

	# grab info about how to find a dataset and read it in to variable "d"
	source(datasets[i])

	# check if data has actually be read in, only then should we fit the models
	if(is.null(d)){
		# print out which dataset WILL NOT be analyzed
		message(paste0("Skipping ",datasetsName))
	}else{
		# check if data has actually be read in, only then should we fit the models
		if(is.null(d)){
			# print out which dataset WILL NOT be analyzed
			message(paste0("Skipping ",datasetsName))
		}else{
			# print out which dataset WILL be analyzed
			message(paste0("Fitting ",datasetsName))

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

			# save a copy of the raw data in case we need it for bootstrapping
			d.orig <- d

			#############################################
			# fit all the functional response models
			# NOTE: optimization is on log-transformed values
			#############################################	 

			# Do data need to be bootstrapped?
			if("Nconsumed1.mean" %in% colnames(d)){
				boot.reps <- 250
			}else{
				boot.reps <- 1
			}

			# fit all model formulations separately
			locals <- list()
			for(modeltype in which.models){
				if(grepl("Hybrid",modeltype)){
					link.funcs <- c("identity", "exp")
				}else{
					link.funcs <- c("identity")
				}

				for(lll in link.funcs){


				message(paste0(" ",modeltype," ",lll),appendLF=FALSE)
				# we will perform all fits boot.reps different times
				local.fits <- foreach(b=1:boot.reps) %dopar% {
					# some fits don't work due to wonkiness in the data so we'll just plow forward when that happens
					bad.fit <- TRUE
			  		while(bad.fit){
						if(any(grepl("[.]mean$",colnames(d.orig)))){
							d <- bootstrap.data(d.orig, this.study$replacement)
						}
			    
			    		if(lll=="identity"){
							# attempt to fit the model and abort if the fit fails for some reason	    	
				    		local.fit <- try(fit.holling.like(d, s=this.study, modeltype=modeltype))
				    	}else{
				    		local.fit <- try(fit.holling.like(d, s=this.study, modeltype=modeltype, phi.transform=exp))
						}

				    	if(!inherits(local.fit, "try-error")){
				    		bad.fit <- FALSE
				    	}
				    }
				    message(".",appendLF=FALSE)
				    # pb$tick(tokens = list(modeltype = modeltype))

				    local.fit
				}

				# pb$finished()
				# we made it out of the loop somewhat miraculously
				message(paste0(" Finished"))

				# create container for the parameter estimates
				local.boots <- make.array(local.fits[[1]], boot.reps)

				# create container for the AIC of the fits
				local.AICs <- lapply(local.fits, AIC)

				# create container for the RMSD of the fits
				local.RMSDs <- lapply(local.fits, resid.metric, metric = 'RMSD')
				
				# create container for the MAD of the fits
				local.MADs <- lapply(local.fits, resid.metric, metric = 'MAD')

				# scrape out the parameter estimates
				for(b in 1:boot.reps){
					local.boots[,,b] <- mytidy(local.fits[[b]])
				}

				# ~~~~~~~~~~~~~~~~~~~~
				# Summarize bootstraps
				# ~~~~~~~~~~~~~~~~~~~~
				local.ests <- apply(local.boots, c(1,2), summarize.boots)

				# pb$terminate()

				# save the key stuff
				locals[[paste(modeltype, lll)]] <- list(
					fit=local.fits[[1]],
					ests=local.ests,
					AICs=local.AICs,
					RMSDs=local.RMSDs,
					MADs=local.MADs
				)
			}
			}

		  	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# save the (last) fits, bootstraps summaries, and some data aspects
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			ffr.fit <- list(
		  		study.info = c(
		  			datasetName = datasetsName,
		  			datadir = datadir,
		  			sample.size = nrow(d),
		  	        # this.study,
		  	        data=d
		  	    ),
				fits = lapply(locals, function(x) x$fit),
				estimates = lapply(locals, function(x) x$ests),
				AICs = lapply(locals, function(x) x$AICs),
				RMSDs = lapply(locals, function(x) x$RMSDs),
				MADs = lapply(locals, function(x) x$MADs)
			)

			# Save the data set fit
			saveRDS(
				ffr.fit,
				file=paste0('../../../results/R/OnePredTwoPrey_fits/', datasetsName,'.Rdata')
			)
		}
	}
}
