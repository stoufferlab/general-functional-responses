
# read in the simplify the raw data
datadir <- 'Reeve_1997'
filename <- 'Tdubius_functional_response.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# rename to standard column names used in fitting code
d <- rawdata[,c("n", "p", "na")]
colnames(d) <- c("Nprey", "Npredator", "Nconsumed")

# DEBUG
# grab units from paper
d$Time <- 1

# WARNING non-integers and NA values muck things up
d <- d[!is.na(d$Nconsumed),]
d$Nconsumed <- round(d$Nconsumed)

# DEBUG initial estimate of attack rate in Holling Type I
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# DEBUG initial estimate of attack rate in Hassell-Varley
# NOTE: optimization is on log-transformed values
x0.rd <- c(-3,0)
