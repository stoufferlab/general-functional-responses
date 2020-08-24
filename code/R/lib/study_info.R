
# convert the textual dataset details to a list of useful variables
study.info <- function(details, datadir, datatype){
	# select only the row of interest
	if(datatype == "One_Predator_One_Prey"){
		details <- subset(
			details,
			Dataset_Folder==datadir & Prey_Richness==1
		)[1,]
	}else{
		details <- subset(
			details,
			Dataset_Folder==datadir & Prey_Richness==2
		)[1,]
	}

	# determine whether or not there are P-1 "predators" interfering
	if(details$Predator_Density_or_Count == "Count"){
		Pminus1 <- TRUE
	}
	else{
		Pminus1 <- FALSE
	}

	# does dataset contain (some) means that need to be bootstraped
	if(details$Original_Means_Compilation == "Original"){
		bootstrap <- FALSE
	}
	else{
		bootstrap <- TRUE
	}

	# determine whether or not the study was used by DeLong and Vasseur
	if(grepl("DeLong",details$PriorUse)){
		delong <- TRUE
	}
	else{
		delong <- FALSE
	}

	# determine whether or not the study has predators or parasites
	if(details$Predator.Parasitoid == "Predator"){
		predator <- TRUE
	}
	else{
		if(details$Predator.Parasitoid == "Parasitoid"){
			predator <- FALSE
		}else{
			predator <- NA
		}
	}

	# define the type of experiment
	if(is.na(details$With_Prey_Replacement)){
		replacement <- NA
	}else{
		if(details$With_Prey_Replacement){
			replacement <- TRUE
		}
		else{
			if(!predator && details$Parasitoid_Type == "Non-discriminating"){
				replacement <- TRUE
			}else{
				replacement <- FALSE
			}
		}
	}

	# scrape the time units so that we can put everything on a common temporal scale
	timeunits <- as.character(details$TimeUnits)

	# scrape the bibtex citation for the study
	cite <- as.character(details$CitationKey)
  
	# scrape the way we obtained the data
	datasource <- as.character(details$DataSource)

	# scrape the bibtex citation for the dataset
	datacitation <- as.character(details$DataCitationKey)

	# scrape the source (Fig or table) of the data (for extracted datasets0
	datafigtablesource <- ifelse(
		as.character(details$DataSource)=='Extracted',
		as.character(details$FigTableSource),
		'-'
	)

	# has this dataset been deemed okay to post to the repo
	ok2post <- as.character(details$PostToRepo)

	# put all the info we need into a list
	rt <- list(
		dataname=datadir,
		Pminus1=Pminus1,
		bootstrap=bootstrap,
		delong=delong,
		predator=predator,
		replacement=replacement,
		timeunits=timeunits,
		cite=cite,
		datasource=datasource,
		datafigtablesource=datafigtablesource,
		datacitation=datacitation,
		ok2post=ok2post
	)
	return(rt)
}
