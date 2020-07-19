rm(list = ls())
# set to FALSE if you want to watch messages in real time
# or TRUE to have them silently saved to file instead.
sinkMessages <- TRUE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
	Sys.getenv("LOGNAME"),
	stouffer = '~/Dropbox/Projects/GenFuncResp/Data',
	marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a few utility functions
source('../lib/study_info.R')
source('../lib/read_data.R')
source('../lib/bootstrap_data.R')
source('../lib/mytidySumm.R')
source('../lib/AA_method.R')
source('../lib/set_params.R')
source('../lib/resid_metrics.R')
source('../lib/plot_coefs.R')
source('../lib/holling_method_one_predator_one_prey.R') # may throw ignorable warning and takes a while to load because of C++ compiling
source('../lib/ratio_method_one_predator_one_prey.R') # may throw ignorable warning and takes a while to load because of C++ compiling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(progress)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# master list of datasets
datasets <- list.files('./Dataset_Code', pattern=".R$", full.names=TRUE, include.dirs=FALSE)

# define the models which are to be fit
holling.like.models <- c(
	"Holling.I",
	"Holling.II",
	"Beddington.DeAngelis",
	"Crowley.Martin",
	"Stouffer.Novak.I"
)

ratio.like.models <- c(
	"Ratio",
	"Arditi.Ginzburg",
	"Hassell.Varley",
	"Arditi.Akcakaya",
	"Arditi.Akcakaya.Method.2"
)

# set the random seed so that bootstrapping is reliable
# generated one integer between 1 and 100000 with Random Integer Generator at random.org
# Timestamp: 2020-07-16 03:52:42 UTC
set.seed(49801)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's start analyzing!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
		# print out which dataset WILL be analyzed
		message(paste0("Fitting ",datasetsName))

		# grab info from the google doc
		this.study <- study.info(datadir)

		# tranform data into terms of hours
		if(!is.null(d$Time)){
			d$Time <- switch(
				this.study$timeunits,
				Seconds = d$Time / 3600.,
				Minutes = d$Time / 60.,
				Hours = d$Time,
				Days = d$Time * 24,
				Unavailable = rep(1, nrow(d))
			)
		}

		#############################################
		# fit all the functional response models
		#############################################

		# start capturing the progress and warning messages
		if(sinkMessages){
			options(warn=1) # provide more than just the base info level
			Mesgs <- file(paste0('../../../results/R/OnePredOnePrey_ErrorLog/', datasetsName, '_ErrorLog.txt'), open='wt')
			sink(Mesgs, type="message")
		}

		# save original data in case we need to bootstrap it
		d.orig <- d

		# Do data need to be bootstrapped?
		if("Nconsumed.mean" %in% colnames(d)){
			boot.reps <- 250
		}else{
			boot.reps <- 1
		}

		# create a progress bar that shows how far along the fitting is
		# the number of reps depends on whether or not we will fit with the AAMethod
		pb <- progress_bar$new(
			format = "  bootstrapping [:bar] :percent eta: :eta",
			total = boot.reps,
			show_after = 0,
			force = TRUE,
			clear = FALSE
		)
		pb$tick(0)

		# temporary storage of for model fits
		bootstrap.fits <- list()
		for(modeltype in c(holling.like.models,ratio.like.models)){
			if(modeltype != "Arditi.Akcakaya.Method.2" || okay4AAmethod(d)){
				bootstrap.fits[[modeltype]] <- list()
			}
		}

		# perform 1 or many bootstrapped fits
		for(b in 1:boot.reps){
			# some bootstrapped data fails for reasons hard to determine
			bad.fit <- TRUE
			while(bad.fit){
				# generate bootstrapped data if necessary
				if("Nconsumed.mean" %in% colnames(d.orig)){
					d <- bootstrap.data(d.orig, this.study$replacement)
				}

				# attempt to fit all models and catch an error if any fit fails
				success <- try({
					for(modeltype in c(holling.like.models,ratio.like.models)){
						if(modeltype != "Arditi.Akcakaya.Method.2" || okay4AAmethod(d)){
							# attempt to fit the model; abort and re-bootstrap if the fit fails
							if(modeltype %in% holling.like.models){
								bootstrap.fits[[modeltype]][[b]] <- fit.holling.like(d, this.study, modeltype)
							}else{
								if(modeltype != "Arditi.Akcakaya.Method.2"){
									bootstrap.fits[[modeltype]][[b]] <- fit.ratio.like(d, this.study, modeltype)
								}else{
									bootstrap.fits[[modeltype]][[b]] <- AAmethod(d, this.study$replacement)
								}
							}
						}
					}
				})

				# the fit succeeded we can continue to the next bootstrap
				if(!inherits(success, "try-error")){
					bad.fit <- FALSE
					pb$tick()
				}
			}
		}

		for(modeltype in c(holling.like.models,ratio.like.models)){
			if(modeltype != "Arditi.Akcakaya.Method.2" || okay4AAmethod(d)){
				# create container for the parameter estimates
				if(modeltype!="Arditi.Akcakaya.Method.2"){
					# ~~~~~~~~~~~~~~~~~~~~
					# Summarize bootstraps
					# ~~~~~~~~~~~~~~~~~~~~
					# scrape out the parameter estimates					
					local.boots <- make.array(bootstrap.fits[[modeltype]][[1]], boot.reps)
					for(b in 1:boot.reps){
						local.boots[,,b] <- mytidy(bootstrap.fits[[modeltype]][[b]])
					}

					# get out the estimates into their own object
					local.ests <- as.array(apply(local.boots, c(1,2), summarize.boots))

					# create container for the logLik of the fits
					local.lls <- summarize.boots(sapply(bootstrap.fits[[modeltype]], logLik))

					# create container for the AIC of the fits
					local.AICs <- summarize.boots(sapply(bootstrap.fits[[modeltype]], AIC))

					# create container for the AICc of the fits
					local.AICcs <- summarize.boots(sapply(bootstrap.fits[[modeltype]], AICc))

					# create container for the RMSD of the fits
					local.RMSDs <- summarize.boots(sapply(bootstrap.fits[[modeltype]], resid.metric, metric = 'RMSD'))

					# create container for the MAD of the fits
					local.MADs <- summarize.boots(sapply(bootstrap.fits[[modeltype]], resid.metric, metric = 'MAD'))

					# save the key stuff (including the first fit)
					bootstrap.fits[[modeltype]] <- list(
						fit=bootstrap.fits[[modeltype]][[1]],
						boots=local.boots,
						ests=local.ests,
						lls=local.lls,
						AICs=local.AICs,
						AICcs=local.AICcs,
						RMSDs=local.RMSDs,
						MADs=local.MADs
					)
				}else{
					# scrape out the parameter estimates
					local.boots <- make.array(bootstrap.fits[[modeltype]][[1]]$estimates, boot.reps)
					for(b in 1:boot.reps){
						local.boots[,,b] <- bootstrap.fits[[modeltype]][[b]]$estimates
					}
					# save the key stuff (which for this model is only a subset of the above)
					bootstrap.fits[[modeltype]] <- list(
						fit=bootstrap.fits[[modeltype]][[1]],
						boots=local.boots,
						ests=local.ests
					)
				}
			}
		}

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# save the (first) fit, bootstraps summaries, and some data aspects
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ffr.fit <- list(
			study.info = c(
				datasetName = datasetsName,
				datadir = datadir,
				sample.size = nrow(d),
				data=d.orig,
				this.study
			),
			fits = lapply(bootstrap.fits, function(x) x$fit),
			boots = lapply(bootstrap.fits, function(x) x$boots),
			estimates = lapply(bootstrap.fits, function(x) x$ests),
			LL = lapply(bootstrap.fits, function(x) x$lls),
			AIC = lapply(bootstrap.fits, function(x) x$AICs),
			AICc = lapply(bootstrap.fits, function(x) x$AICcs),
			RMSD = lapply(bootstrap.fits, function(x) x$RMSDs),
			MAD = lapply(bootstrap.fits, function(x) x$MADs)
		)
		
		# Save the data set fit monster object
		saveRDS(
			ffr.fit,
			file=paste0('../../../results/R/OnePredOnePrey_fits/', datasetsName,'.Rdata')
		)

		# close open streams, etc
		if(sinkMessages){
			sink(type="message")
			close(Mesgs)
			options(warn=0)
			readLines(paste0('../../../results/R/OnePredOnePrey_ErrorLog/', datasetsName, '_ErrorLog.txt'))
		}
	}
}
