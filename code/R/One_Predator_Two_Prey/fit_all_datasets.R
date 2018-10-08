rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
	Sys.getenv("LOGNAME"),
	stouffer = '../../../dropbox_data/Data',
	MARKTOUPDATETHIS = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a few utility functions
source('../lib/study_info.R')
source('../lib/bootstrap_data.R')
# source('../LogLikelihoods/AA_method.R')
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
# datasets <- c("./Dataset_Code/vucetich_2002_isleroyale_whole2014.R")

# fit everything on a dataset by dataset basis
for(i in 1:length(datasets)){
	# loads the data into data frame 'd' and specifies data-specific parameters
	source(datasets[i])

	# grab some info from the google doc
	this.study <- study.info(datadir)

	if(!grepl("T", this.study$runswith)){
		message(paste0("No to ",datasets[i]))
	}else{
		# print out which dataset is being analyzed
		message(paste0("Yes to ",datasets[i]))

		#############################################
		# fit all the functional response models
		# NOTE: optimization is on log-transformed values
		#############################################	 

		# Do data need to be bootstrapped? If so, save the raw data frame d as d.org.
		d.orig <- d
		if("Nconsumed1.mean" %in% colnames(d)){
			boot.reps <- 2
			d <- bootstrap.data(d.orig, this.study$response)
		}else{
			boot.reps <- 1
		}

		# for initial fits to produce dimensions for containers
		ffr.hollingI <- fit.holling.like(d, s=this.study, modeltype="Holling I")
		ffr.hollingII.specialist.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Specialist")
		# ffr.hollingII.specialist.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Generalist")
		# ffr.hollingII.generalist.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Specialist")
		ffr.hollingII.generalist.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Generalist")
		# ffr.hollingII.specialist.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Hybrid")
		# ffr.hollingII.generalist.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Hybrid")
		# ffr.hollingII.hybrid.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Specialist")
		# ffr.hollingII.hybrid.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Generalist")
		ffr.hollingII.hybrid.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Hybrid")

		# # Containers for estimates
		boots.HT.I <- make.array(ffr.hollingI, boot.reps)
		boots.HT.II.SS <- make.array(ffr.hollingII.specialist.specialist, boot.reps)
		boots.HT.II.GG <- make.array(ffr.hollingII.generalist.generalist, boot.reps)
		boots.HT.II.HH <- make.array(ffr.hollingII.hybrid.hybrid, boot.reps)
		# boots.BD <- make.array(ffr.bd, boot.reps)
		# boots.CM <- make.array(ffr.cm, boot.reps)
		# boots.SN.I <- make.array(ffr.sn1, boot.reps)
		# # boots.SN.Numer <- make.array(ffr.sn2, boot.reps)
		# # boots.SN.III <- make.array(ffr.sn3, boot.reps)
		# # boots.HV <- make.array(ffr.hv, boot.reps)
		# # boots.AG <- make.array(ffr.ag, boot.reps)
		# # boots.AA <- make.array(ffr.aa, boot.reps)
		# # if(okay4AAmethod(d)){
		# # 	boots.AA2 <- array(NA, dim=c(dim(fit.AAmethod$estimates), boot.reps))
		# #     dimnames(boots.AA2)[c(1,2)] <- dimnames(fit.AAmethod$estimates)
	 # 	# }
		  
	  	for(b in 1:boot.reps){
			if(any(grepl("[.]mean$",colnames(d.orig)))){
				d <- bootstrap.data(d.orig, this.study$response)
			}
	    
			ffr.hollingI <- fit.holling.like(d, s=this.study, modeltype="Holling I")
			ffr.hollingII.specialist.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Specialist")
			# ffr.hollingII.specialist.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Generalist")
			# ffr.hollingII.generalist.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Specialist")
			ffr.hollingII.generalist.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Generalist")
			# ffr.hollingII.specialist.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Specialist Hybrid")
			# ffr.hollingII.generalist.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Generalist Hybrid")
			# ffr.hollingII.hybrid.specialist <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Specialist")
			# ffr.hollingII.hybrid.generalist <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Generalist")
			ffr.hollingII.hybrid.hybrid <- fit.holling.like(d, s=this.study, modeltype="Holling II Hybrid Hybrid")

	 #    	# source('fit_holling_like_nobounds.R')
	 #    	fit.hollingI <- fit.holling.like(d, this.study, "Holling I")
		# 	ffr.hollingII <- fit.holling.like(d, this.study, "Holling II")
		# 	ffr.bd <- fit.holling.like(d, this.study, "Beddington-DeAngelis")
		# 	ffr.cm <- fit.holling.like(d, this.study, "Crowley-Martin")
		# 	ffr.sn1 <- fit.holling.like(d, this.study, "Stouffer-Novak I")
	 #    	# ifelse(okay4AAmethod(d), fit.AAmethod <- AAmethod(d,expttype), fit.AAmethod <- NA)
	    	
			boots.HT.I[,,b] <- mytidy(ffr.hollingI)
		 	boots.HT.II.SS[,,b] <- mytidy(ffr.hollingII.specialist.specialist)
		 	boots.HT.II.GG[,,b] <- mytidy(ffr.hollingII.generalist.generalist)
		 	boots.HT.II.HH[,,b] <- mytidy(ffr.hollingII.hybrid.hybrid)
		# 	boots.BD[,,b] <- mytidy(ffr.bd)
		# 	boots.CM[,,b] <- mytidy(ffr.cm)
		# 	boots.SN.I[,,b] <- mytidy(ffr.sn1)
		# 	# boots.SN.Numer[,,b] <- mytidy(ffr.sn2)
		# 	# boots.SN.III[,,b] <- mytidy(ffr.sn3)
		# 	# boots.HV[,,b] <- mytidy(ffr.hv)
		# 	# boots.AG[,,b] <- mytidy(ffr.ag)
		# 	# boots.AA[,,b] <- mytidy(ffr.aa)
		# 	# if(okay4AAmethod(d)){ boots.AA2[,,b] <- fit.AAmethod$estimates  }
		  
			# # ~~~~~~~~~~~~~~~~~~~~
			# # Summarize bootstraps
			# # ~~~~~~~~~~~~~~~~~~~~
			HT.I.ests <- as.array(apply(boots.HT.I, c(1,2), summarize.boots))
			HT.II.SS.ests <- as.array(apply(boots.HT.II.SS, c(1,2), summarize.boots))
			HT.II.GG.ests <- as.array(apply(boots.HT.II.GG, c(1,2), summarize.boots))
			HT.II.HH.ests <- as.array(apply(boots.HT.II.HH, c(1,2), summarize.boots))
			# HT.II.ests <- as.array(apply(boots.HT.II, c(1,2), summarize.boots))
			# BD.ests <- as.array(apply(boots.BD, c(1,2), summarize.boots))
			# CM.ests <- as.array(apply(boots.CM, c(1,2), summarize.boots))
			# SN.I.ests <- as.array(apply(boots.SN.I, c(1,2), summarize.boots))
			# # SN.Numer.ests <- as.array(apply(boots.SN.Numer,c(1,2), summarize.boots))
			# # # SN.III.ests <- as.array(apply(boots.SN.III,c(1,2), summarize.boots))
			# # HV.ests <- as.array(apply(boots.HV,c(1,2), summarize.boots))	  
			# # AG.ests <- as.array(apply(boots.AG,c(1,2), summarize.boots))
			# # AA.ests <- as.array(apply(boots.AA,c(1,2), summarize.boots))
			# # if(okay4AAmethod(d)){ 
			# #     AA2.ests <- as.array(apply(boots.AA2,c(1,2), summarize.boots))
			# # }else{
			# # 	AA2.ests <- NA
			# # }
		}
  
	 #  	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# # save the (last) fits and some data aspects
		# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		ffr.fits[[datasets[i]]] <- list(
	  		study.info = c(
	  			datadir = datadir,
	  			sample.size = nrow(d),
	  	        # this.study,
	  	        data=d
	  	    ),
			fits = c(
				Holling.Type.II.Specialist.Specialist = ffr.hollingII.specialist.specialist,
				# Holling.Type.II.Specialist.Generalist = ffr.hollingII.specialist.generalist,
				# Holling.Type.II.Generalist.Specialist = ffr.hollingII.generalist.specialist,
				Holling.Type.II.Generalist.Generalist = ffr.hollingII.generalist.generalist,
				# Holling.Type.II.Specialist.Hybrid = ffr.hollingII.specialist.hybrid,
				# Holling.Type.II.Generalist.Hybrid = ffr.hollingII.generalist.hybrid,
				# Holling.Type.II.Hybrid.Specialist = ffr.hollingII.hybrid.specialist,
				# Holling.Type.II.Hybrid.Generalist = ffr.hollingII.hybrid.generalist,
				Holling.Type.II.Hybrid.Hybrid = ffr.hollingII.hybrid.hybrid,
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

		break
	}

	# source('plot_phi_denom.R')
	# plot.AAmethod(fit.AAmethod)
	# break
}

# # save the mega container which includes all FR fits
# save(ffr.fits,file='../../../results/R/ffr.fits_OnePredTwoPrey.Rdata')

# # generate a quick and dirty plot of the phi_denom parameters of the SNI model
# source('plot_phi_denom.R')
