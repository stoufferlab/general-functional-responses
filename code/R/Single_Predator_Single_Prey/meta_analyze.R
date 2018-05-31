
# specify where the data files are located
dropboxdir <- '../../dropbox_data/Data'

# a few utility functions
source('study_info.R')
source('bootstrap_data.R')

# cobble together a master list of things to analyze (yes, this is very clunky right now)
datasets <- do.call(rbind,list(
	# c("Ganjisaffar 2017","ganjisaffar_2017.R"),
	# c("Nilsson 2004","nilsson_2004.R"),
	c("Katz 1985","katz_1985.R"), # WARNING this datasets gives errors because of the bootstrap
	c("Chan 2017 Coyotes Hare","chan_2017_coyotes_hare.R"),
	c("Chan 2017 Coyotes Squirrels","chan_2017_coyotes_squirrels.R"),
	c("Chan 2017 Lynx Hare","chan_2017_lynx_hare.R"),
	c("Chan 2017 Lynx Squirrels","chan_2017_lynx_squirrels.R"),
	c("Chant 1966","chant_1966.R"),
	c("Chong 2006","chong_2006.R"),
	c("Crowley 1989","crowley_1989.R"),
	c("Edwards 1961 Nasonia Musca","edwards_1961_Nasonia-Musca.R"),
	c("Edwards 1961 Trichogramma Sitotroga 1","edwards_1961_Trichogramma-Sitotroga_1.R"),
	c("Edwards 1961 Trichogramma Sitotroga 2","edwards_1961_Trichogramma-Sitotroga_2.R"),
	c("Eveleigh 1982 Adeg adult","eveleigh_1982_Adeg_adult.R"),
	c("Eveleigh 1982 Adeg deut","eveleigh_1982_Adeg_deut.R"),
	c("Eveleigh 1982 Adeg proto","eveleigh_1982_Adeg_proto.R"),
	c("Eveleigh 1982 Pper adult","eveleigh_1982_Pper_adult.R"),
	c("Eveleigh 1982 Pper deut","eveleigh_1982_Pper_deut.R"),
	c("Eveleigh 1982 Pper proto","eveleigh_1982_Pper_proto.R"),
	c("Hassan 1976 Aglom","hassan_1976_Aglom.R"),
	c("Hassan 1976 Breg","hassan_1976_Breg.R"),
	c("Hassan 1976 Ppup","hassan_1976_Ppup.R"),
	c("Hossie 2016 Clumped","hossie_2016_clumped.R"),
	c("Hossie 2016 Even","hossie_2016_even.R"),
	c("Kratina 2009","kratina_2009.R"),
	c("Lang 2012 Poe 10C","lang_2012_Poe_10C.R"),
	c("Lang 2012 Poe 20C","lang_2012_Poe_20C.R"),
	c("Lang 2012 Pter 10C","lang_2012_Pter_10C.R"),
	c("Lang 2012 Pter 20C","lang_2012_Pter_20C.R"),
	c("Long 2012a","long_2012a.R"),
	c("Mertz 1968","mertz_1968.R"),
	c("Omkar 2004","omkar_2004.R"),
	c("Puscak 2018","pusack_2018.R"),
	c("Reeve 1997","reeve_1997.R"),
	c("Salt 1974","salt_1974.R"),
	c("Walde 1984","walde_1984.R"),
	c("Wasserman 2016 Bluegill","wasserman_2016_bluegill.R"),
	c("Wasserman 2016 Mouthbrooder","wasserman_2016_mouthbrooder.R"),
	c("Wasserman 2016 Tilapia","wasserman_2016_tilapia.R")
))

# for all of the above go through and fit all holling-like and ratio-dependent-like FRs
ffr.fits <- list()
for(i in 1:nrow(datasets)){
	# message(datasets[i,1])

	# loads the data and specifies all parameters
	source(paste0("Dataset_Code/",datasets[i,2]))

	# DEBUG should arguably do bootstrapping here

	# fits all the functional response models
	source('fit_holling_like.R')
	source('fit_ratio_dependent.R')

	# save the fits and some data aspects to a "convenient" list
	ffr.fits[[datasets[i,1]]] <- list(
		"Holling Type I" = ffr.hollingI,
		"Holling Type II" = ffr.hollingII,
		"Beddington-DeAngelis" = ffr.bd,
		"Crowley-Martin" = ffr.cm,
		"Stouffer-Novak I" = ffr.sn1,
		# "Stouffer-Novak Numer" = ffr.sn2,
		"Stouffer-Novak II" = ffr.sn3,
		"Hassell-Varley" = ffr.hv,
		"Arditi-Ginzburg" = ffr.ag,
		"Arditi-Akcakaya" = ffr.aa,
		"sample.size" = nrow(d),
		"expttype" = expttype,
		"Pminus1" = Pminus1,
		"datadir" = datadir,
		"data" = d
	)
	# break
}

##########################################
# some heinous plotting code appears below
##########################################

studies <- names(ffr.fits)

# sort by sample size?
if(FALSE){
	studies <- names(sort(sapply(ffr.fits, function(x) x[["sample.size"]])))
}

# prepare two plots, one of phi and one of m
par(mfrow=c(1,2))

# plot of phi (denom)
plot(
	1:length(ffr.fits) ~ unlist(lapply(ffr.fits, function(x){max(-4,min(4,coef(x[["Stouffer-Novak I"]])["phi_denom"]))})),
	type='n',
	ylab="Dataset",
	xlab="Phi [ Crowley-Martin = -4 ; Beddington-DeAngelis=4 ]",
	xlim=c(-4,4),
	yaxt='n'
)
axis(side=2, at=1:length(ffr.fits), labels=FALSE)

i <- 1
for(nn in studies){
	x <- ffr.fits[[nn]]
	mm <- coef(x[["Stouffer-Novak I"]])["phi_denom"]
	if(mm < 4 & mm > -4){
		se <- coef(summary(x[["Stouffer-Novak I"]]))["phi_denom","Std. Error"]
		if(!is.nan(se)){
			segments((mm-se), i, (mm+se), i)
		}
	}else{
		mm <- max(-4,min(4,mm))
	}
	points(y=c(i),x=c(mm))
	i <- i + 1
}

# plot of m (ratio-dependent exponent)
if(TRUE){
	plot(
		1:length(ffr.fits) ~ unlist(lapply(ffr.fits, function(x){exp(coef(x[["Arditi-Akcakaya"]])["exponent"])})),
		type='n',
		ylab="",
		xlab="Ratio-dependence m",
		yaxt='n',
		xlim=c(0,4),
	)
	axis(side=2, at=1:length(ffr.fits), labels=FALSE)

	i <- 1
	for(nn in studies){
		x <- ffr.fits[[nn]]
		mm <- coef(x[["Arditi-Akcakaya"]])["exponent"]
		se <- coef(summary(x[["Arditi-Akcakaya"]]))["exponent","Std. Error"]
		if(!is.nan(se)){
			segments(exp(mm-se), i, exp(mm+se), i)
		}
		points(y=c(i),x=c(exp(mm)))
		i <- i + 1
	}

	abline(v=1)
}

# leftovers are hit tv series
# print(AICtab(ffr.hollingI, ffr.hollingII, ffr.bd, ffr.cm, ffr.sn1, ffr.sn2, ffr.sn3, weights=TRUE))
