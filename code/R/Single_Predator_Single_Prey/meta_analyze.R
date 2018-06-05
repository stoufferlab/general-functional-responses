
# specify where the data files are located
dropboxdir <- '../../../dropbox_data/Data'

# a few utility functions
source('study_info.R')
source('bootstrap_data.R')

# cobble together a master list of things to analyze (yes, this is very clunky right now)
datasets <- list.files('./Dataset_Code',full.names=TRUE)
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)

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

	# should arguably do multiple bootstraps here
	if("Nconsumed.mean" %in% colnames(d)){
		d <- bootstrap.data(d, expttype)
	}

	#############################################
	# fit all the functional response models
	# should turn these into functions for sure
	#############################################

	# DEBUG initial estimate of attack rate in Holling Type I
	# NOTE: optimization is on log-transformed values
	x0.hl <- c(attack=log(0.01))
	source('fit_holling_like.R')

	# DEBUG initial estimate of attack rate in Hassell-Varley
	# NOTE: optimization is on log-transformed values
	x0.rd <- c(attack=log(0.01),exponent=log(1))
	source('fit_ratio_dependent.R')

	# save the fits and some data aspects to a "convenient" list
	ffr.fits[[datasets[i]]] <- list(
		"data" = d,
		"sample.size" = nrow(d),
		"expttype" = expttype,
		"Pminus1" = Pminus1,
		"datadir" = datadir,
		"study.info" = this.study,
		"Holling Type I" = ffr.hollingI,
		"Holling Type II" = ffr.hollingII,
		"Beddington-DeAngelis" = ffr.bd,
		"Crowley-Martin" = ffr.cm,
		"Stouffer-Novak I" = ffr.sn1,
		# "Stouffer-Novak Numer" = ffr.sn2,
		"Stouffer-Novak II" = ffr.sn3,
		"Hassell-Varley" = ffr.hv,
		"Arditi-Ginzburg" = ffr.ag,
		"Arditi-Akcakaya" = ffr.aa
	)
	# break
}

##########################################
# some heinous plotting code appears below
##########################################

source('plot_coefs.R')

# sort by sample size?
if(TRUE){
	studies <- names(sort(sapply(ffr.fits, function(x) x[["sample.size"]])))
	ffr.fits <- ffr.fits[studies]
}

ffr.fits.delong <- ffr.fits[unlist(lapply(ffr.fits, function(x) x$study.info$delong))]

# prepare some plots of the fits
# DEBUG work on the confidence intervals aspect
par(mfrow=c(2,2))

# plot.coefs(ffr.fits,
# 	"Holling Type II",
# 	"attack",
# 	plot.SEs=TRUE,
# 	ilink=exp,
# 	xlab="Holling Type II Attack Rate",
# 	ylab="Dataset",
# 	xlim=c(0,3)
# )

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
