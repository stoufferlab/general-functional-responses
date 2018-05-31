
study.info <- function(dataname){
	# grab some info from the google doc
	require(RCurl)
	googledoc<-read.csv(
		text=getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vQcykYqM8Pkmgrlmp9S2jorZZEOlZ14a0AINRuDc2Y_29f6dTR9ojhOOBV2rcattJO5LXA5ATVn_nK6/pub?gid=0&single=true&output=csv"),
		header=T,
		skip=1,
		sep=",",
		na.strings = c("*", "NA")
	)
	googledoc <- googledoc[googledoc$Dataset_Folder == datadir & !is.na(googledoc$Dataset_Folder),]

	# define the type of experiment
	if(googledoc[1,"With_Prey_Replacement"] | is.na(googledoc[1,"With_Prey_Replacement"])){
		expttype <- "replacement"
	}
	else{
		expttype <- "integrated"
	}

	# determine whether or not there are P-1 "predators" interfering
	if(googledoc[1,"Predator_Density_or_Count"] == "Count"){
		Pminus1 <- TRUE
	}
	else{
		Pminus1 <- FALSE
	}

	# define the type of experiment
	if(googledoc[1,"Count_per_Predator"] | is.na(googledoc[1,"Count_per_Predator"])){
		perpredator <- TRUE
	}
	else{
		perpredator <- FALSE
	}

	rt <- list(
		expttype=expttype,
		Pminus1=Pminus1
	)
}
