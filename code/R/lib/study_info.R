
study.info <- function(dataname){
  
	# grab some info from the google doc
	# require(RCurl)
	# masterlist<-read.csv(
#     text=getURL(
#       "https://docs.google.com/spreadsheets/d/e/2PACX-1vQcykYqM8Pkmgrlmp9S2jorZZEOlZ14a0AINRuDc2Y_29f6dTR9ojhOOBV2rcattJO5LXA5ATVn_nK6/pub?gid=0&single=true&output=csv",
#       .opts = list(followlocation = TRUE)
#     ),
#     header=T,
#     skip=1,
#     sep=",",
# 		na.strings = c("*", "NA")
# 	)
	
	masterlist <- read.csv(paste(dropboxdir,'GenFunResp-PotentialData.csv',sep="/"),
	                       header=T,
	                       skip=1,
	                       sep=",",
	                       na.strings = c("*", "NA"))
	
	masterlist <- masterlist[masterlist$Dataset_Folder == dataname & !is.na(masterlist$Dataset_Folder),]
	
	# determine whether or not there are P-1 "predators" interfering
	if(masterlist[1,"Predator_Density_or_Count"] == "Count"){
		Pminus1 <- TRUE
	}
	else{
		Pminus1 <- FALSE
	}

	# does dataset contain (some) means that need to be bootstraped
	if(masterlist[1,"Original_Means_Compilation"] == "Original"){
		bootstrap <- FALSE
	}
	else{
		bootstrap <- TRUE
	}

	# determine whether or not the study was used by DeLong and Vasseur
	if(grepl("DeLong",masterlist[1,"PriorUse"])){
		delong <- TRUE
	}
	else{
		delong <- FALSE
	}

	# determine whether or not the study has predators or parasites
	if(masterlist[1,"Predator.Parasitoid"] == "Predator"){
		predator <- TRUE
	}
	else{
		if(masterlist[1,"Predator.Parasitoid"] == "Parasitoid"){
			predator <- FALSE
		}else{
			predator <- NA
		}
	}

	# define the type of experiment
	if(is.na(masterlist[1,"With_Prey_Replacement"])){
		replacement <- NA
	}else{
		if(masterlist[1,"With_Prey_Replacement"]){
			replacement <- TRUE
		}
		else{
			if(!predator && masterlist[1,"Parasitoid_Type"] == "Non-discriminating"){
				replacement <- TRUE
			}else{
				replacement <- FALSE
			}
		}
	}

	# determine whether or not the data can be fit by the holling-like and/or ratio-dependent codes
	runswith <- as.character(masterlist[1,"RunsWith"])

	# scrape the time units so that we can put everything on a common temporal scale
	timeunits <- as.character(masterlist[1,"TimeUnits"])
  
  # scrape the bibtex citation for the study
  cite <- as.character(masterlist[1,"CitationKey"])
  
	# put all the info we need into a list
	rt <- list(
	  dataname=dataname,
		Pminus1=Pminus1,
		bootstrap=bootstrap,
		delong=delong,
		predator=predator,
		replacement=replacement,
		runswith=runswith,
		timeunits=timeunits,
		cite=cite
	)
}
