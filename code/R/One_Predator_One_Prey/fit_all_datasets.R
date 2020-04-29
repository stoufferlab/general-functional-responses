rm(list = ls())
# set to FALSE if you want to watch messages in real time
# or TRUE to have them silently saved to file instead.
sinkMessages <- TRUE
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
source('../lib/set_params.R')
source('../lib/RMSD.R')
source('../lib/plot_coefs.R')
source('../lib/holling_method_one_predator_one_prey.R') # takes a while to load because of C++ compiling
source('../lib/ratio_method_one_predator_one_prey.R') # takes a while to load because of C++ compiling
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

# check to see which datasets have been fit (and thus on't bother refitting them)
# datasets.fitted <- list.files('../../../results/R/OnePredOnePrey_fits/', full.names=FALSE, include.dirs=FALSE)
# datasets.fitted <- gsub('*data$','',datasets.fitted)
# datasets <- datasets[!gsub('./Dataset_Code/','',datasets)%in%datasets.fitted]

# select focal dataset for testing
# datasets <- c("./Dataset_Code/Walde_1984.R")  # Occasional Hessian problem
# datasets <- datasets[1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's start analyzing!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit everything on a dataset by dataset basis
for(i in 1:length(datasets)){

	# loads the data into data frame 'd' and specifies data-specific parameters
	source(datasets[i])

	# grab info from the google doc
	this.study <- study.info(datadir)

	datasetsName <- sub('*.R$','', sub('*./Dataset_Code/','', datasets[i]))

	#############################################
	# fit all the functional response models
	# NOTE: optimization is on log-transformed values for some of the parameters
	#############################################

	# H: holling-like, R: ratio-like, T: test set (or combinations thereof)
	# if(!grepl("R", this.study$runswith)){
	if(!grepl("H|R", this.study$runswith)){
		# print out which dataset is not being analyzed
		message(paste0("No to ",datasetsName))
	}else{
		# print out which dataset is being analyzed
		message(paste0("Yes to ",datasetsName))

		# start capturing the progress and warning messages
		if(sinkMessages){
			options(warn=1) # provide more than just the base info level
			Mesgs <- file(paste0('../../../results/R/OnePredOnePrey_ErrorLog/', datasetsName, '_ErrorLog.txt'), open='wt')
			sink(Mesgs, type="message")
		}

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

		# save original data in case we need to bootstrap it
		d.orig <- d

		# Do data need to be bootstrapped?
		if("Nconsumed.mean" %in% colnames(d)){
			boot.reps <- 250
		}else{
			boot.reps <- 1
		}

		# create a progress bar that shows how far along the bootstrapping is
		pb <- progress_bar$new(
			format = "  bootstrapping [:bar] :percent eta: :eta",
			total = boot.reps*10,
			show_after = 0,
			force = TRUE,
			clear = FALSE
		)
		pb$tick(0)

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
					ffr.hollingI <- fit.holling.like(d, this.study, "Holling.I")
					pb$tick()
					ffr.hollingII <- fit.holling.like(d, this.study, "Holling.II")
					pb$tick()
					ffr.bd <- fit.holling.like(d, this.study, "Beddington.DeAngelis")
					pb$tick()
					ffr.cm <- fit.holling.like(d, this.study, "Crowley.Martin")
					pb$tick()
					ffr.sn1 <- fit.holling.like(d, this.study, "Stouffer.Novak.I")
					pb$tick()
				}else{
					for(tick in 1:5) pb$tick()
				}

				ffr.ratio <- ffr.ag <- ffr.hv <- ffr.aa <- ffr.aam <- array(NA,c(1,1))
				if(grepl("R", this.study$runswith)){
					ffr.ratio <- fit.ratio.like(d, this.study, "Ratio")
					pb$tick()
					ffr.ag <- fit.ratio.like(d, this.study, "Arditi.Ginzburg")
					pb$tick()
					ffr.hv <- fit.ratio.like(d, this.study, "Hassell.Varley")
					pb$tick()
					ffr.aa <- fit.ratio.like(d, this.study, "Arditi.Akcakaya")
					pb$tick()
					if(okay4AAmethod(d)){
						ffr.aam <- AAmethod(d,this.study$replacement)
					}
					pb$tick()
				}else{
					for(tick in 1:5) pb$tick()
				}
			})

			if(inherits(success, "try-error")){
				pb$update((b-1)/boot.reps)
			} else {
				if(b == 1){
					# create containers for log likelihood values of all fits
					ll.hollingI <- ll.hollingII <- ll.bd <- ll.cm <- ll.sn1 <- rep(NA,boot.reps)
					# ll.sn2 <- ll.sn3 <- vector("list",boot.reps)
					ll.ratio <- ll.ag <- ll.hv <- ll.aa <- rep(NA,boot.reps)

					# create containers for AIC values of all fits
					aic.hollingI <- aic.hollingII <- aic.bd <- aic.cm <- aic.sn1 <- rep(NA,boot.reps)
					# aic.sn2 <- aic.sn3 <- vector("list",boot.reps)
					aic.ratio <- aic.ag <- aic.hv <- aic.aa <- rep(NA,boot.reps)

					# create containers for AICc values of all fits
					aicc.hollingI <- aicc.hollingII <- aicc.bd <- aicc.cm <- aicc.sn1 <- rep(NA,boot.reps)
					# aicc.sn2 <- aicc.sn3 <- vector("list",boot.reps)
					aicc.ratio <- aicc.ag <- aicc.hv <- aicc.aa <- rep(NA,boot.reps)

					# create containers for RMSD values of all fits
					rmsd.hollingI <- rmsd.hollingII <- rmsd.bd <- rmsd.cm <- rmsd.sn1 <- rep(NA,boot.reps)
					# RMSD.sn2 <- RMSD.sn3 <- vector("list",boot.reps)
					rmsd.ratio <- rmsd.ag <- rmsd.hv <- rmsd.aa <- rep(NA,boot.reps)

					# create containers for parameter estimates using first bootstrap as a template
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

				# Save the AICc values for the rep
				if(grepl("H", this.study$runswith)){
					ll.hollingI[[b]] <- logLik(ffr.hollingI)
					ll.hollingII[[b]] <- logLik(ffr.hollingII)
					ll.bd[[b]] <- logLik(ffr.bd)
					ll.cm[[b]] <- logLik(ffr.cm)
					ll.sn1[[b]] <- logLik(ffr.sn1)
					# ll.sn2[[b]] <- logLik(ffr.sn2)
					# ll.sn3[[b]] <- logLik(ffr.sn3)

					aic.hollingI[[b]] <- AIC(ffr.hollingI)
					aic.hollingII[[b]] <- AIC(ffr.hollingII)
					aic.bd[[b]] <- AIC(ffr.bd)
					aic.cm[[b]] <- AIC(ffr.cm)
					aic.sn1[[b]] <- AIC(ffr.sn1)
					# aic.sn2[[b]] <- AIC(ffr.sn2)
					# aic.sn3[[b]] <- AIC(ffr.sn3)

					aicc.hollingI[[b]] <- AICc(ffr.hollingI)
					aicc.hollingII[[b]] <- AICc(ffr.hollingII)
					aicc.bd[[b]] <- AICc(ffr.bd)
					aicc.cm[[b]] <- AICc(ffr.cm)
					aicc.sn1[[b]] <- AICc(ffr.sn1)
					# aicc.sn2[[b]] <- AICc(ffr.sn2)
					# aicc.sn3[[b]] <- AICc(ffr.sn3)

					rmsd.hollingI[[b]] <- RMSD(ffr.hollingI)
					rmsd.hollingII[[b]] <- RMSD(ffr.hollingII)
					rmsd.bd[[b]] <- RMSD(ffr.bd)
					rmsd.cm[[b]] <- RMSD(ffr.cm)
					rmsd.sn1[[b]] <- RMSD(ffr.sn1)
					# RMSD.sn2[[b]] <- RMSD(d, ffr.sn2, this.study, 'Stouffer.Novak.II')
					# RMSD.sn3[[b]] <- RMSD(d, ffr.sn3, this.study, 'Stouffer.Novak.III')
				}

				if(grepl("R", this.study$runswith)){
					ll.ratio[[b]] <- logLik(ffr.ratio)
					ll.ag[[b]] <- logLik(ffr.ag)
					ll.hv[[b]] <- logLik(ffr.hv)
					ll.aa[[b]] <- logLik(ffr.aa)

					aic.ratio[[b]] <- AIC(ffr.ratio)
					aic.ag[[b]] <- AIC(ffr.ag)
					aic.hv[[b]] <- AIC(ffr.hv)
					aic.aa[[b]] <- AIC(ffr.aa)

					aicc.ratio[[b]] <- AICc(ffr.ratio)
					aicc.ag[[b]] <- AICc(ffr.ag)
					aicc.hv[[b]] <- AICc(ffr.hv)
					aicc.aa[[b]] <- AICc(ffr.aa)

					rmsd.ratio[[b]] <- RMSD(ffr.ratio)
					rmsd.ag[[b]] <- RMSD(ffr.ag)
					rmsd.hv[[b]] <- RMSD(ffr.hv)
					rmsd.aa[[b]] <- RMSD(ffr.aa)
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

				# update the bootstrap reps counter
				b <- b + 1
			}
		}

		# ~~~~~~~~~~~~~~~~~~~~
		# Summarize bootstraps
		# ~~~~~~~~~~~~~~~~~~~~
		ests.hollingI <- ests.hollingII <- ests.bd <- ests.cm <- ests.sn1 <- NA
		LL.hollingI <- LL.hollingII <- LL.bd <- LL.cm <- LL.sn1 <- NA
		AIC.hollingI <- AIC.hollingII <- AIC.bd <- AIC.cm <- AIC.sn1 <- NA
		AICc.hollingI <- AICc.hollingII <- AICc.bd <- AICc.cm <- AICc.sn1 <- NA
		RMSD.hollingI <- RMSD.hollingII <- RMSD.bd <- RMSD.cm <- RMSD.sn1 <- NA
		if(grepl("H", this.study$runswith)){
			ests.hollingI <- as.array(apply(boots.hollingI, c(1,2), summarize.boots))
			ests.hollingII <- as.array(apply(boots.hollingII, c(1,2), summarize.boots))
			ests.bd <- as.array(apply(boots.bd, c(1,2), summarize.boots))
			ests.cm <- as.array(apply(boots.cm, c(1,2), summarize.boots))
			ests.sn1 <- as.array(apply(boots.sn1, c(1,2), summarize.boots))
			# ests.sn2 <- as.array(apply(boots.sn2,c(1,2), summarize.boots))
			# ests.sn3 <- as.array(apply(boots.sn3,c(1,2), summarize.boots))

			LL.hollingI <- summarize.boots(ll.hollingI)
			LL.hollingII <- summarize.boots(ll.hollingII)
			LL.bd <- summarize.boots(ll.bd)
			LL.cm <- summarize.boots(ll.cm)
			LL.sn1 <- summarize.boots(ll.sn1)
			# LL.sn2 <- summarize.boots(ll.sn2)
			# LL.sn3 <- summarize.boots(ll.sn3)

			AIC.hollingI <- summarize.boots(aic.hollingI)
			AIC.hollingII <- summarize.boots(aic.hollingII)
			AIC.bd <- summarize.boots(aic.bd)
			AIC.cm <- summarize.boots(aic.cm)
			AIC.sn1 <- summarize.boots(aic.sn1)
			# AIC.sn2 <- summarize.boots(aic.sn2)
			# AIC.sn3 <- summarize.boots(aic.sn3)

			AICc.hollingI <- summarize.boots(aicc.hollingI)
			AICc.hollingII <- summarize.boots(aicc.hollingII)
			AICc.bd <- summarize.boots(aicc.bd)
			AICc.cm <- summarize.boots(aicc.cm)
			AICc.sn1 <- summarize.boots(aicc.sn1)
			# AICc.sn2 <- summarize.boots(aicc.sn2)
			# AICc.sn3 <- summarize.boots(aicc.sn3)

			RMSD.hollingI <- summarize.boots(rmsd.hollingI)
			RMSD.hollingII <- summarize.boots(rmsd.hollingII)
			RMSD.bd <- summarize.boots(rmsd.bd)
			RMSD.cm <- summarize.boots(rmsd.cm)
			RMSD.sn1 <- summarize.boots(rmsd.sn1)
			# RMSD.sn2 <- summarize.boots(RMSD.sn2)
			# RMSD.sn3 <- summarize.boots(RMSD.sn3)
		}

		ests.ratio <- ests.ag <- ests.hv <- ests.aa <- ests.aam <- NA
		LL.ratio <- LL.ag <- LL.hv <- LL.aa <- LL.aam <- NA
		AIC.ratio <- AIC.ag <- AIC.hv <- AIC.aa <- AIC.aam <- NA
		AICc.ratio <- AICc.ag <- AICc.hv <- AICc.aa <- AICc.aam <- NA
		RMSD.ratio <- RMSD.ag <- RMSD.hv <- RMSD.aa <- RMSD.aam <- NA
		if(grepl("R", this.study$runswith)){
			ests.ratio <- as.array(apply(boots.ratio,c(1,2), summarize.boots))
			ests.ag <- as.array(apply(boots.ag,c(1,2), summarize.boots))
			ests.hv <- as.array(apply(boots.hv,c(1,2), summarize.boots))
			ests.aa <- as.array(apply(boots.aa,c(1,2), summarize.boots))
			if(okay4AAmethod(d)){
				ests.aam <- as.array(apply(boots.aam,c(1,2), summarize.boots))
			}

			LL.ratio <- summarize.boots(ll.ratio)
			LL.ag <- summarize.boots(ll.ag)
			LL.hv <- summarize.boots(ll.hv)
			LL.aa <- summarize.boots(ll.aa)

			AIC.ratio <- summarize.boots(aic.ratio)
			AIC.ag <- summarize.boots(aic.ag)
			AIC.hv <- summarize.boots(aic.hv)
			AIC.aa <- summarize.boots(aic.aa)

			AICc.ratio <- summarize.boots(aicc.ratio)
			AICc.ag <- summarize.boots(aicc.ag)
			AICc.hv <- summarize.boots(aicc.hv)
			AICc.aa <- summarize.boots(aicc.aa)

			RMSD.ratio <- summarize.boots(rmsd.ratio)
			RMSD.ag <- summarize.boots(rmsd.ag)
			RMSD.hv <- summarize.boots(rmsd.hv)
			RMSD.aa <- summarize.boots(rmsd.aa)
		}

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# save study info, the (last) fits, bootstrap estimates and estimate summaries
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# ffr.fits[[datasets[i]]] <- list(
		ffr.fit <- list(
			study.info = c(
				datasetName = datasetsName,
				datadir = datadir,
				sample.size = nrow(d),
				data=d.orig,
				this.study
			),
			fits = list(
				Holling.I = ffr.hollingI,
				Holling.II = ffr.hollingII,
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
				Holling.I = boots.hollingI,
				Holling.II = boots.hollingII,
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
			LL = list(
				Holling.I = LL.hollingI,
				Holling.II = LL.hollingII,
				Beddington.DeAngelis = LL.bd,
				Crowley.Martin = LL.cm,
				Stouffer.Novak.I = LL.sn1,
				# Stouffer.Novak.II = LL.sn2,
				# Stouffer.Novak.III = LL.sn3,
				Ratio = LL.ratio,
				Arditi.Ginzburg = LL.ag,
				Hassell.Varley = LL.hv,
				Arditi.Akcakaya = LL.aa
			),
			AIC = list(
				Holling.I = AIC.hollingI,
				Holling.II = AIC.hollingII,
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
			AICc = list(
				Holling.I = AICc.hollingI,
				Holling.II = AICc.hollingII,
				Beddington.DeAngelis = AICc.bd,
				Crowley.Martin = AICc.cm,
				Stouffer.Novak.I = AICc.sn1,
				# Stouffer.Novak.II = AICc.sn2,
				# Stouffer.Novak.III = AICc.sn3,
				Ratio = AICc.ratio,
				Arditi.Ginzburg = AICc.ag,
				Hassell.Varley = AICc.hv,
				Arditi.Akcakaya = AICc.aa
			),
			RMSD = list(
				Holling.I = RMSD.hollingI,
				Holling.II = RMSD.hollingII,
				Beddington.DeAngelis = RMSD.bd,
				Crowley.Martin = RMSD.cm,
				Stouffer.Novak.I = RMSD.sn1,
				# Stouffer.Novak.II = RMSD.sn2,
				# Stouffer.Novak.III = RMSD.sn3,
				Ratio = RMSD.ratio,
				Arditi.Ginzburg = RMSD.ag,
				Hassell.Varley = RMSD.hv,
				Arditi.Akcakaya = RMSD.aa
			),
			estimates = list(
				Holling.I = ests.hollingI,
				Holling.II = ests.hollingII,
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
		
		# Save the data set fit
		saveRDS(
			ffr.fit,
			file=paste0('../../../results/R/OnePredOnePrey_fits/', datasetsName,'.Rdata')
		)

		if(sinkMessages){
			sink(type="message")
			close(Mesgs)
			options(warn=0)
			readLines(paste0('../../../results/R/OnePredOnePrey_ErrorLog/', datasetsName, '_ErrorLog.txt'))
		}
	}
}
# sink(type="message")
# close(Mesgs)

# save a mega container
ffr.fits <- bundle_fits('../../../results/R/OnePredOnePrey_fits')
save(
	ffr.fits,
	file='../../../results/R/OnePredOnePrey_ffr.fits.Rdata'
)

