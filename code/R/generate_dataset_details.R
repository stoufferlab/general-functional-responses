
require(RCurl)
# grab the information from the google doc
masterlist<-read.csv(
	text=getURL(
		"https://docs.google.com/spreadsheets/d/e/2PACX-1vQcykYqM8Pkmgrlmp9S2jorZZEOlZ14a0AINRuDc2Y_29f6dTR9ojhOOBV2rcattJO5LXA5ATVn_nK6/pub?gid=0&single=true&output=csv",
		.opts = list(followlocation = TRUE)
	),
	header=T,
	skip=1,
	sep=",",
	na.strings = c("*", "NA")
	)

# only keep datasets that we used
masterlist <- masterlist[which(masterlist[,1]),]

# there only a subset of columns we need and use
tokeep <- c(
	"Dataset_Folder",
	"Variable_Predator_Abundance",
	"Prey_Richness",
	"Predator_Density_or_Count",
	"Original_Means_Compilation",
	"PriorUse",
	"Predator.Parasitoid",
	"With_Prey_Replacement",
	"Parasitoid_Type",
	"TimeUnits",
	"CitationKey",
	"DataSource",
	"DataCitationKey",
	"FigTableSource",
	"PostToRepo"
)

# create a reduced data frame
strippedlist <- masterlist[,tokeep]
rownames(strippedlist) <- 1:nrow(strippedlist)

# write out a csv to be used from within the repo
write.table(strippedlist,file="dataset_details.csv")

