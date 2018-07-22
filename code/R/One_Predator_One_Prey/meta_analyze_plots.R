##########################################
# some heinous plotting code appears below
##########################################
load(file='../../../results/R/ffr.fits_OnePredOnePrey.Rdata')


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
