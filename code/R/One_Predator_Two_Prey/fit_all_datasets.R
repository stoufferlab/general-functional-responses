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
registerDoParallel(cores=7)
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

# # DEBUG: for testing only
# datasets <- c("./Dataset_Code/Iyer_1996_Br.R","./Dataset_Code/Iyer_1996_Bp.R","./Dataset_Code/Iyer_1996_Bc.R")

# fit everything on a dataset by dataset basis
for(i in 1:length(datasets)){
	# loads the data into data frame 'd' and specifies data-specific parameters
	source(datasets[i])

	# save a copy of the raw data in case we need it for bootstrapping
	d.orig <- d

	# grab some info from the google doc
	this.study <- study.info(datadir)

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

		library(progress)
		pb <- progress_bar$new(
			format = "  bootstrapping [:bar] :percent eta: :eta",
			total = boot.reps,
			show_after = 0,
			force = TRUE,
			clear = FALSE
		)

		b <- 1
	  	while(b <= boot.reps){
			if(any(grepl("[.]mean$",colnames(d.orig)))){
				d <- bootstrap.data(d.orig, this.study$response)
			}
	    
	    	# DEBUG: sometimes the fits fail; should we allow the code to skip these and keep on keepin on in that case?

	    	# fit a suite of functional response models
	    	success <- try({
				ffr.hollingI <- fit.holling.like(d, s=this.study, modeltype="Holling I")
				ffr.hollingII.SS <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Specialist")
				# ffr.hollingII.specialist.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Generalist")
				# ffr.hollingII.generalist.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Specialist")
				ffr.hollingII.GG <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Generalist")
				# ffr.hollingII.specialist.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Hybrid")
				# ffr.hollingII.generalist.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Hybrid")
				# ffr.hollingII.hybrid.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Specialist")
				# ffr.hollingII.hybrid.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Generalist")
				ffr.hollingII.HH <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Hybrid")
			})

	    	if(!inherits(success, "try-error")){
				# create containers for the parameter estimates
				if(b == 1){
					boots.HT.I <- make.array(ffr.hollingI, boot.reps)
					boots.HT.II.SS <- make.array(ffr.hollingII.SS, boot.reps)
					boots.HT.II.GG <- make.array(ffr.hollingII.GG, boot.reps)
					boots.HT.II.HH <- make.array(ffr.hollingII.HH, boot.reps)
				}

				# add fits to the containers
				boots.HT.I[,,b] <- mytidy(ffr.hollingI)
			 	boots.HT.II.SS[,,b] <- mytidy(ffr.hollingII.SS)
			 	boots.HT.II.GG[,,b] <- mytidy(ffr.hollingII.GG)
			 	boots.HT.II.HH[,,b] <- mytidy(ffr.hollingII.HH)

			 	pb$tick()
			 	b <- b + 1
			 }
		 }

		# ~~~~~~~~~~~~~~~~~~~~
		# Summarize bootstraps
		# ~~~~~~~~~~~~~~~~~~~~
		HT.I.ests <- as.array(apply(boots.HT.I, c(1,2), summarize.boots))
		HT.II.SS.ests <- as.array(apply(boots.HT.II.SS, c(1,2), summarize.boots))
		HT.II.GG.ests <- as.array(apply(boots.HT.II.GG, c(1,2), summarize.boots))
		HT.II.HH.ests <- as.array(apply(boots.HT.II.HH, c(1,2), summarize.boots))
  
	  	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# save the (last) fits and some data aspects
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ffr.fits[[datasets[i]]] <- list(
	  		study.info = c(
	  			datadir = datadir,
	  			sample.size = nrow(d),
	  	        # this.study,
	  	        data=d
	  	    ),
			fits = c(
				Holling.Type.II.Specialist.Specialist = ffr.hollingII.SS,
				# Holling.Type.II.Specialist.Generalist = ffr.hollingII.specialist.generalist,
				# Holling.Type.II.Generalist.Specialist = ffr.hollingII.generalist.specialist,
				Holling.Type.II.Generalist.Generalist = ffr.hollingII.GG,
				# Holling.Type.II.Specialist.Hybrid = ffr.hollingII.specialist.hybrid,
				# Holling.Type.II.Generalist.Hybrid = ffr.hollingII.generalist.hybrid,
				# Holling.Type.II.Hybrid.Specialist = ffr.hollingII.hybrid.specialist,
				# Holling.Type.II.Hybrid.Generalist = ffr.hollingII.hybrid.generalist,
				Holling.Type.II.Hybrid.Hybrid = ffr.hollingII.HH,
				Holling.Type.I = ffr.hollingI
	        ),
			estimates = list(
			    Holling.Type.II.Specialist.Specialist = HT.II.SS.ests,
				# Holling.Type.II.Specialist.Generalist = ffr.hollingII.specialist.generalist,
				# Holling.Type.II.Generalist.Specialist = ffr.hollingII.generalist.specialist,
				Holling.Type.II.Generalist.Generalist = HT.II.GG.ests,
				# Holling.Type.II.Specialist.Hybrid = ffr.hollingII.specialist.hybrid,
				# Holling.Type.II.Generalist.Hybrid = ffr.hollingII.generalist.hybrid,
				# Holling.Type.II.Hybrid.Specialist = ffr.hollingII.hybrid.specialist,
				# Holling.Type.II.Hybrid.Generalist = ffr.hollingII.hybrid.generalist,
				Holling.Type.II.Hybrid.Hybrid = HT.II.HH.ests,
				Holling.Type.I = HT.I.ests
			)
		)
	}

	# source('plot_phi_denom.R')
	# plot.AAmethod(fit.AAmethod)
	# break
}

# save the mega container which includes all FR fits
save(ffr.fits,file='../../../results/R/ffr.fits_OnePredTwoPrey.Rdata')

# # generate a quick and dirty plot of the phi_denom parameters of the SNI model
# source('plot_phi_denom.R')
