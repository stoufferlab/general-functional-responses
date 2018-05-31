
# read in the simplify the raw data
datadir <- 'Pusack_2018'
filename <- 'Pusack et al. 2018 Data for M Novak.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# rename to standard column names used in fitting code
d <- rawdata[,c("drill.abundance", "oyster.abundance", "total.no.oysters.consumed", "days")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")

# DEBUG initial estimate of attack rate in Holling Type I
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# DEBUG initial estimate of attack rate in Hassell-Varley
# NOTE: optimization is on log-transformed values
x0.rd <- c(-3,0)

