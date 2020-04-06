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
source('../lib/holling_method_one_predator_two_prey.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################
library(doParallel)
registerDoParallel(cores=6)
####################################

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
	#DEBUG
	out <- rbind(out, c(AIC(fit),0,0,0))
	rownames(out)[nrow(out)] <- "AIC"
	#ENDDEBUG
	out
}

make.array <- function(ffr.fit,boot.reps){
	t.ffr.fit <- mytidy(ffr.fit)
	out <- array(NA, dim=c(dim(t.ffr.fit), boot.reps))
	dimnames(out) <- dimnames(t.ffr.fit)
	out
}

summarize.boots <- function(x){
	c(mean=mean(x, na.rm=TRUE), quantile(x,c(0.025,0.16,0.5,0.84,0.975), na.rm=TRUE), n=sum(!is.na(x)))
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

# # # DEBUG: for testing only
# datasets <- datasets[[6]] #c("./Dataset_Code/Kalinkat_2011_Anch.R") #,"./Dataset_Code/zzz_Buckel_2000_small.R")

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
			Unavailable = rep(1, nrow(d))
		)
	}

	# save a copy of the raw data in case we need it for bootstrapping
	d.orig <- d

	if(!grepl("H", this.study$runswith)){
		message(paste0("No to ",datasets[i]))
	}else{
		# print out which dataset is being analyzed
		message(paste0("Yes to ",datasets[i]))

		#############################################
		# fit all the functional response models
		# NOTE: optimization is on log-transformed values
		#############################################	 

		# Do data need to be bootstrapped?
		if("Nconsumed1.mean" %in% colnames(d)){
			boot.reps <- 100
		}else{
			boot.reps <- 1
		}

		# library(progress)

		# fit all model formulations separately
		locals <- list()
		for(modeltype in which.models){
			# pb <- progress::progress_bar$new(
			# 	format = "  :modeltype [:bar] :percent eta: :eta",
			# 	total = boot.reps,
			# 	show_after = 0,
			# 	force = TRUE,
			# 	clear = FALSE
			# )
			# pb$tick(0, tokens = list(modeltype = modeltype))

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
			locals[[paste(modeltype, lll)]] <- list(fit=local.fits[[1]], ests=local.ests)
		}
		}

	  	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# save the (last) fits, bootstraps summaries, and some data aspects
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ffr.fits[[datasets[i]]] <- list(
	  		study.info = c(
	  			datadir = datadir,
	  			sample.size = nrow(d),
	  	        # this.study,
	  	        data=d
	  	    ),
			fits = lapply(locals, function(x) x$fit),
			estimates = lapply(locals, function(x) x$ests)
		)

	}

	# source('plot_phi_denom.R')
	# plot.AAmethod(fit.AAmethod)
	# break
}

# save the mega container which includes all FR fits
save(ffr.fits,file='../../../results/R/OnePredTwoPrey_ffr.fits.Rdata')

# # generate a quick and dirty plot of the phi_denom parameters of the SNI model
# source('plot_phi_denom.R')
