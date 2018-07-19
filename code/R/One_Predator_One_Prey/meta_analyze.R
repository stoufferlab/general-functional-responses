
# specify where the data files are located
dropboxdir <- '../../../dropbox_data/Data' # Stouffer
dropboxdir <- '~/Dropbox/Research/Projects/GenFuncResp/Data' # Novak


# a few utility functions
source('study_info.R')
source('bootstrap_data.R')
source('../LogLikelihoods/AA_method.R')
# library(devtools)
# devtools::install_github("bbolker/broom")
library(broom) # for summarizing bootstrap fits

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cobble together a master list of things to analyze (yes, this is very clunky right now)
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)
datasets <- grep("zzz", datasets, invert=TRUE, value=TRUE)

# for all of the above go through and fit all holling-like and ratio-dependent-like FRs
ffr.fits <- list()
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
	if("Nconsumed.mean" %in% colnames(d)){
	  d.orig <- d
	  boot.reps <- 3
	} else{boot.reps <- 1}
	
	#############################################
	# fit all the functional response models
	# should turn these into functions for sure
	#############################################	  
	ffr.fits[[datasets[i]]] <- vector('list', boot.reps)
	
	for(b in 1:boot.reps){
    if(b > 1){
	    d <- bootstrap.data(d.orig, expttype)
    }
  	# DEBUG initial estimate of attack rate in Holling Type I
  	# NOTE: optimization is on log-transformed values
  	x0.hl <- c(attack=log(0.01))
  
  	# try({
  	source('fit_holling_like_nobounds.R')
  	ifelse(okay4AAmethod(d), fit.AAmethod <- AAmethod(d,expttype), fit.AAmethod <- NA)
  
  	# # DEBUG initial estimate of attack rate in Hassell-Varley
  	# # NOTE: optimization is on log-transformed values
  	# x0.rd <- c(attack=log(0.01),exponent=log(1))
  	# source('fit_ratio_dependent.R')
  
  	# save the fits and some data aspects to a "convenient" list

  	ffr.fits[[datasets[i]]][[b]] <- list(
    	study.info =    c(data=d,
    	                     this.study,
    	                    sample.size = nrow(d),
    	                    Pminus1 = Pminus1,
    	                    datadir = datadir),
  		Holling.Type.I = list(fit=ffr.hollingI, estimates=tidy(ffr.hollingI)),
  		Holling.Type.II = list(fit=ffr.hollingI, estimates=tidy(ffr.hollingI)),
  		# Beddington.DeAngelis = list(fit=ffr.bd, estimates=tidy(ffr.bd)),
  		# Crowley.Martin = list(fit=ffr.cm, estimates=tidy(ffr.cm)),
  		Stouffer.Novak.I = list(fit=ffr.sn1, estimates=tidy(ffr.sn1)),
  		# Stouffer.Novak.Numer = list(fit=ffr.sn2, estimates=tidy(ffr.sn2)),
  		# Stouffer.Novak.III = list(fit=ffr.sn3, estimates=tidy(ffr.sn3)),
  		# Hassell.Varley = list(fit=ffr.hv, estimates=tidy(ffr.hv)),
  		# Arditi.Ginzburg = list(fit=ffr.ag, estimates=tidy(ffr.ag)),
  		# Arditi.Akcakaya = list(fit=ffr.aa, estimates=tidy(ffr.aa)),
  		Arditi.Akcakaya.Method.2 = fit.AAmethod
  	)
  
  	# })
  
  	# source('plot_phi_denom.R')
  	# plot.AAmethod(fit.AAmethod)
  	# break
	}
	}

# ~~~~~~~~~~~~~~~~~~~~
# Summarize bootstraps
# ~~~~~~~~~~~~~~~~~~~~
for(i in 1:length(datasets)){
  fits <- ffr.fits[[datasets[i]]]
  fitsb <- fits[[1]]
}

# test <- list(A = list(a = c(a = 10, `2` = 20, `3` = 30, `4` = 72),
#                        b = c(`1` = 15, `2` = 9, `3` = 7)),
#              B = list(a = c(A = 11, B = 12, C = 13),
#                        b = c(X = 14, Y = 15, Z = 16)))
# sapply(test, function(x) x[['b']])

x<-data.frame(t(sapply(fits, function(x) x$Holling.Type.I$estimates)))













##########################################
# some heinous plotting code appears below
##########################################

source('plot_coefs.R') # Function is missing from Git

# sort by sample size?
if(TRUE){
	studies <- names(sort(sapply(ffr.fits, function(x) x[["sample.size"]])))
	ffr.fits <- ffr.fits[studies]
}

ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]

# prepare some plots of the fits
# DEBUG work on the confidence intervals aspect
# par(mfrow=c(2,2))

plot.coefs(ffr.fits,
	"Stouffer-Novak I",
	"phi_denom",
	plot.SEs=FALSE,
	ilink=identity,
	xlab="Effect of foraging on interfering",
	ylab="Dataset",
	xlim=c(-4,4)
)

# # plot.coefs(ffr.fits,
# # 	"Beddington-DeAngelis"
# # 	"interference",
# # 	plot.SEs=TRUE,
# # 	ilink=exp,
# # 	xlab="Beddington-DeAngelis Interference",
# # 	# ylab="Dataset"
# # )

# plot.coefs(ffr.fits,
# 	"Hassell-Varley",
# 	"exponent",
# 	plot.SEs=TRUE,
# 	ilink=exp,
# 	xlab="Hassell-Varley Exponent",
# 	ylab="",
# 	xlim=c(0,4)
# )

plot.coefs(ffr.fits,
	"Arditi-Akcakaya",
	"handling",
	plot.SEs=TRUE,
	ilink=exp,
	xlab="Arditi-Akcakaya Handling Time",
	ylab="Dataset",
	xlim=c(0,1)
)

plot.coefs(ffr.fits,
	"Arditi-Akcakaya",
	"exponent",
	plot.SEs=TRUE,
	ilink=exp,
	xlab="Arditi-Akcakaya Exponent",
	ylab=""#,
	# xlim=c(0,4)
)


plot.coefs(ffr.fits.delong,
	"Arditi-Akcakaya",
	"handling",
	plot.SEs=TRUE,
	ilink=exp,
	xlab="Arditi-Akcakaya Handling Time",
	ylab="Dataset",
	xlim=c(0,1)
)

plot.coefs(ffr.fits.delong,
	"Arditi-Akcakaya",
	"exponent",
	plot.SEs=TRUE,
	ilink=exp,
	xlab="Arditi-Akcakaya Exponent",
	ylab=""#,
	# xlim=c(0,4)
)

# leftovers are hit tv series
# print(AICtab(ffr.hollingI, ffr.hollingII, ffr.bd, ffr.cm, ffr.sn1, ffr.sn2, ffr.sn3, weights=TRUE))
