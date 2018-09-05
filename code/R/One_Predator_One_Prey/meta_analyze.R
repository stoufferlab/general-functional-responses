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
source('study_info.R')
source('bootstrap_data.R')
source('../LogLikelihoods/AA_method.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(devtools)
# devtools::install_github("bbolker/broom")
library(broom) # for tidy()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  c(mean=mean(x), quantile(x,c(0.025,0.975)))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# master list of datasets
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)

# remove template files which don't actually read data
datasets <- grep("^template",datasets,invert=TRUE,value=TRUE)

# remove zzz files which are placeholders while a dataset is being cleaned/incorporated
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)

# create a container for the things that get fit
ffr.fits <- vector('list', length(datasets))
names(ffr.fits) <- datasets

for(i in 1:length(datasets)){
	message(datasets[i])

	# loads the data into data frame 'd' and specifies data-specific parameters
	source(datasets[i])

	# grab some info from the google doc
	this.study <- study.info(datadir)
	expttype <- this.study$expttype
	Pminus1 <- this.study$Pminus1

	# # for my debugging
	# print(d$Time)

	# Do data need to be bootstrapped?
	d.orig <- d
	if("Nconsumed.mean" %in% colnames(d)){
	  boot.reps <- 3
	  d <- bootstrap.data(d.orig, expttype)
	} else{boot.reps <- 1}
	
	#############################################
	# fit all the functional response models
	# NOTE: optimization is on log-transformed values
	#############################################	 
	
  # Initial estimate of attack rate in Holling Type I
  # x0.hl <- c(attack=log(0.01))
	NP <- d$Npredator*d$Nprey
	x0.hl <- c(attack=log(coef(lm(d$Nconsumed~0+NP))))
  
  # try({
	# for initial fits to produce dimensions for containers
  source('fit_holling_like_nobounds.R')
  ifelse(okay4AAmethod(d), fit.AAmethod <- AAmethod(d,expttype), fit.AAmethod <- NA)

	# Containers for estimates
	  boots.HT.I <- make.array(ffr.hollingI, boot.reps)
	  boots.HT.II <- make.array(ffr.hollingII, boot.reps)
	  # boots.BD <- make.array(ffr.bd, boot.reps)
	  # boots.CM <- make.array(ffr.cm, boot.reps)
	  boots.SN.I <- make.array(ffr.sn1, boot.reps)
	  # boots.SN.Numer <- make.array(ffr.sn2, boot.reps)
	  # boots.SN.III <- make.array(ffr.sn3, boot.reps)
	  # boots.HV <- make.array(ffr.hv, boot.reps)
	  # boots.AG <- make.array(ffr.ag, boot.reps)
	  # boots.AA <- make.array(ffr.aa, boot.reps)
	  if(okay4AAmethod(d)){
	    boots.AA2 <- array(NA, dim=c(dim(fit.AAmethod$estimates), boot.reps))
	    dimnames(boots.AA2)[c(1,2)] <- dimnames(fit.AAmethod$estimates)
    }
	  
  	for(b in 1:boot.reps){
  	  if("Nconsumed.mean" %in% colnames(d.orig)){
  	  d <- bootstrap.data(d.orig, expttype)
  	  }
    	x0.hl <- c(attack=log(0.01))
    
    	source('fit_holling_like_nobounds.R')
    	ifelse(okay4AAmethod(d), fit.AAmethod <- AAmethod(d,expttype), fit.AAmethod <- NA)
    	
  	  boots.HT.I[,,b] <- mytidy(ffr.hollingI)
  	  boots.HT.II[,,b] <- mytidy(ffr.hollingII)
  	  # boots.BD[,,b] <- mytidy(ffr.bd)
  	  # boots.CM[,,b] <- mytidy(ffr.cm)
  	  boots.SN.I[,,b] <- mytidy(ffr.sn1)
  	  # boots.SN.Numer[,,b] <- mytidy(ffr.sn2)
  	  # boots.SN.III[,,b] <- mytidy(ffr.sn3)
  	  # boots.HV[,,b] <- mytidy(ffr.hv)
  	  # boots.AG[,,b] <- mytidy(ffr.ag)
  	  # boots.AA[,,b] <- mytidy(ffr.aa)
  	  if(okay4AAmethod(d)){ boots.AA2[,,b] <- fit.AAmethod$estimates  }
  	}
	  
	  # ~~~~~~~~~~~~~~~~~~~~
	  # Summarize bootstraps
	  # ~~~~~~~~~~~~~~~~~~~~
	  HT.I.ests <- as.array(apply(boots.HT.I,c(1,2), summarize.boots))
	  HT.II.ests <- as.array(apply(boots.HT.II,c(1,2), summarize.boots))
# 	  BD.ests <- as.array(apply(boots.BD,c(1,2), summarize.boots))
#   	CM.ests <- as.array(apply(boots.CM,c(1,2), summarize.boots))
	  SN.I.ests <- as.array(apply(boots.SN.I,c(1,2), summarize.boots))
# 	  SN.Numer.ests <- as.array(apply(boots.SN.Numer,c(1,2), summarize.boots))
# 	  SN.III.ests <- as.array(apply(boots.SN.III,c(1,2), summarize.boots))
# 	  HV.ests <- as.array(apply(boots.HV,c(1,2), summarize.boots))	  
# 	  AG.ests <- as.array(apply(boots.AG,c(1,2), summarize.boots))
# 	  AA.ests <- as.array(apply(boots.AA,c(1,2), summarize.boots))
	  if(okay4AAmethod(d)){ 
	    AA2.ests <- as.array(apply(boots.AA2,c(1,2), summarize.boots))
	 } else{AA2.ests <- NA}
	  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# save the (last) fits and some data aspects
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	ffr.fits[[datasets[i]]] <- list(
  		study.info = c(
  			data=d,
  	        this.study,
  	        sample.size = nrow(d),
  	        Pminus1 = Pminus1,
  	        datadir = datadir
  	    ),
		fits = c(
			Holling.Type.I = ffr.hollingI,
		    Holling.Type.II = ffr.hollingI,
          	# Beddington.DeAngelis = ffr.bd,
          	# Crowley.Martin = ffr.cm,
          	Stouffer.Novak.I = ffr.sn1
          	# Stouffer.Novak.Numer = ffr.sn2,
          	# Stouffer.Novak.III = ffr.sn3,
          	# Hassell.Varley = ffr.hv,
          	# Arditi.Ginzburg = ffr.ag,
          	# Arditi.Akcakaya = ffr.aa,
          	# Arditi.Akcakaya.Method.2 = fit.AAmethod
        ),
		estimates = list(
		    Holling.Type.I = HT.I.ests,
		    Holling.Type.II = HT.II.ests,
		    # Beddington.DeAngelis = BD.ests,
		    # Crowley.Martin = CM.ests,
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

	# source('plot_phi_denom.R')
	# plot.AAmethod(fit.AAmethod)
	# break
}

save(ffr.fits,file='../../../results/R/ffr.fits_OnePredOnePrey.Rdata')




