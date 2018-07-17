
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
	googledoc <- googledoc[googledoc$Dataset_Folder == dataname & !is.na(googledoc$Dataset_Folder),]

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

	# does dataset contain (some) means that need to be bootstraped
	if(googledoc[1,"Original_Means_Compilation"] == "Original"){
		bootstrap <- FALSE
	}
	else{
		bootstrap <- TRUE
	}

	# determine whether or not the study was used by DeLong and Vasseur
	if(grepl("DeLong",googledoc[1,"PriorUse"])){
		delong <- TRUE
	}
	else{
		delong <- FALSE
	}

	# determine whether or not the study has predators or parasites
	if(googledoc[1,"Predator.Parasitoid"] == "Predator"){
		predator <- TRUE
	}
	else{
		if(googledoc[1,"Predator.Parasitoid"] == "Parasitoid"){
			predator <- FALSE
		}else{
			predator <- NA
		}
	}

	rt <- list(
		expttype=expttype,
		Pminus1=Pminus1,
		bootstrap=bootstrap,
		delong=delong,
		predator=predator
	)
}
