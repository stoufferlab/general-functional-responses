rm(list = ls())
# set to FALSE if you want to watch messages in real time
# or TRUE to have them silently saved to file instead.
sinkMessages <- TRUE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a few utility functions
source('../lib/bootstrap_data.R')
source('../lib/mytidySumm.R')
source('../lib/plot_coefs.R')
source('../lib/read_data.R')
source('../lib/resid_metrics.R')
source('../lib/study_info.R')
source('../lib/holling_method_one_predator_two_prey.R') # may throw ignorable warning and takes a while to load because of C++ compiling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################
library(doParallel)
registerDoParallel(cores=6)
####################################

# read in the table of dataset details
dataset_details <- read.csv(
	'../../../data/dataset_details.csv'
)

# master list of datasets
datasets <- list.files('./Dataset_Code', pattern=".R$", full.names=TRUE, include.dirs=FALSE)

# define the models which are to be fit
which.models <- c(
	"Holling I",
	"Holling II Specialist Specialist",
	"Holling II Specialist Generalist",
	"Holling II Generalist Specialist",
	"Holling II Generalist Generalist",
	"Holling II Hybrid Hybrid"
	# "Holling II Specialist Hybrid",
	# "Holling II Hybrid Specialist",
	# "Holling II Generalist Hybrid",
	# "Holling II Hybrid Generalist",
)

# set the random seed so that bootstrapping is reliable
# generated one integer between 1 and 100000 with Random Integer Generator at random.org
# Timestamp: 2020-07-16 03:53:52 UTC
set.seed(12114)

# fit everything on a dataset by dataset basis
for(i in seq_along(datasets)){
	# create a short nickname for the dataset
	datasetsName <- sub('*.R$','', sub('*./Dataset_Code/','', datasets[i]))

	# grab info about how to find a dataset
	source(datasets[i])

	# use the above information to read the data into a variable d
	d <- read.data(datadir, filename, "One_Predator_Two_Prey", columns)

	# check if data has actually be read in, only then should we fit the models
	if(is.null(d)){
		# print out which dataset WILL NOT be analyzed
		message(paste0("Skipping ",datasetsName))
	}else{
		# check if data has actually be read in, only then should we fit the models
		if(is.null(d)){
			# print out which dataset WILL NOT be analyzed
			cat(paste0("Skipping ",datasetsName,"\n"))
		}else{
			# print out which dataset WILL be analyzed
			cat(paste0("Fitting ",datasetsName,"\n"))

			# grab info about experimental design, etc
			this.study <- study.info(
				dataset_details,
				datadir,
				"One_Predator_Two_Prey"
			)

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

			# start capturing the progress and warning messages
			if(sinkMessages){
				options(warn=1) # provide more than just the base info level
				Mesgs <- file(paste0('../../../results/R/OnePredTwoPrey_ErrorLog/', datasetsName, '_ErrorLog.txt'), open='wt')
				sink(Mesgs, type="message")
			}

			# Do data need to be bootstrapped?
			if("Nconsumed1.mean" %in% colnames(d)){
				boot.reps <- 250
			}else{
				boot.reps <- 1
			}

			# create a progress bar that shows how far along the fitting is
			pb <- txtProgressBar(
				min=0,
				max=boot.reps
			)

			# fit 1 to many bootstrapped datasets
			bootstrap.fits <- foreach(b=1:boot.reps) %dopar% {
				# some fits don't work due to wonkiness in the data so we'll just plow forward when that happens
				bad.fit <- TRUE
				while(bad.fit){
					if(any(grepl("[.]mean$",colnames(d.orig)))){
						d <- bootstrap.data(d.orig, this.study$replacement)
					}

					# fit all model formulations separately
					local.fits <- list()
					success <- try({
						for(modeltype in which.models){
							if(grepl("Hybrid",modeltype)){
								link.funcs <- c("identity", "exp")
							}else{
								link.funcs <- c("identity")
							}

							for(lll in link.funcs){
								if(lll=="identity"){
									# attempt to fit the model and abort if the fit fails for some reason
									local.fits[[paste(modeltype, lll)]] <- fit.holling.like(d, s=this.study, modeltype=modeltype)
								}else{
									local.fits[[paste(modeltype, lll)]] <- fit.holling.like(d, s=this.study, modeltype=modeltype, phi.transform=exp)
								}
							}
						}
					})
					if(!inherits(success, "try-error")){
						bad.fit <- FALSE
						setTxtProgressBar(pb, b)
					}
				}
				local.fits
			}

			# we made it out of the loop somewhat miraculously
			close(pb)
			# message(paste0(" Finished"))

			# bootstrap fits is organized by bootstrapped data
			# reorganize to be based on models
			locals <- list()
			for(modeltype in names(bootstrap.fits[[1]])){
				model.fits <- list()
				for(b in 1:boot.reps){
					model.fits[[b]] <- bootstrap.fits[[b]][[modeltype]]
				}

				# create container for the parameter estimates
				local.boots <- make.array(model.fits[[1]], boot.reps)

				# create container for the AIC of the fits
				local.AICs <- lapply(model.fits, AIC)

				# create container for the RMSD of the fits
				local.RMSDs <- lapply(model.fits, resid.metric, metric = 'RMSD')
				
				# create container for the MAD of the fits
				local.MADs <- lapply(model.fits, resid.metric, metric = 'MAD')

				# scrape out the parameter estimates
				for(b in 1:boot.reps){
					local.boots[,,b] <- mytidy(model.fits[[b]])
				}

				# ~~~~~~~~~~~~~~~~~~~~
				# Summarize bootstraps
				# ~~~~~~~~~~~~~~~~~~~~
				local.ests <- apply(local.boots, c(1,2), summarize.boots)

				# save the key stuff
				locals[[modeltype]] <- list(
					fit=model.fits[[1]],
					ests=local.ests,
					AICs=local.AICs,
					RMSDs=local.RMSDs,
					MADs=local.MADs
				)
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

			# close open streams, etc
			if(sinkMessages){
				sink(type="message")
				close(Mesgs)
				options(warn=0)
				readLines(paste0('../../../results/R/OnePredTwoPrey_ErrorLog/', datasetsName, '_ErrorLog.txt'))
			}
		}
	}
}
