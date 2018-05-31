
# read in the simplify the raw data
datadir <- 'Wasserman_2016'
filename <- 'Data_for_repository.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# tidy up some column names (was more useful when still using stan)
colnames(rawdata) <- gsub("[.]","_",colnames(rawdata))

# determine which predators were present
predators <- unique(grep(" [+] ",rawdata$Predator,invert=TRUE,value=TRUE))
rawdata[,predators] <- 0

# fill in predator abundances
for(r in 1:nrow(rawdata)){
	prds <- unlist(strsplit(as.character(rawdata$Predator[r])," [+] "))
	for(s in 1:length(prds)){
		rawdata[r,prds[s]] <- rawdata[r,prds[s]] + 1
	}
}

# eliminate the text field with predator data
rawdata$Predator <- NULL

# DEBUG this is what makes the different Wasserman files different
d <- subset(rawdata, bluegill>0 & mouthbrooder==0 & tilapia==0)
d$Npredator <- d$bluegill

# remove extraneous columns and give them the common name format
d <- d[,c("Npredator","supplied_prey","eaten_prey")]
colnames(d) <- c("Npredator","Nprey","Nconsumed")

# DEBUG
# grab units from paper
d$Time <- 1

# DEBUG initial estimate of attack rate in Holling Type I
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# DEBUG initial estimate of attack rate in Hassell-Varley
# NOTE: optimization is on log-transformed values
x0.rd <- c(-3,0)
