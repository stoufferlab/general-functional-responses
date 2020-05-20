
# we might need this here
library(bbmle)

# read in the different fits
ffr.fits <- readRDS(
	file='../../../../results/R/OnePredTwoPrey_ffr.fits.Rdata'
)

# scrape out the dataset names
labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))

# scrape out the AIC values for the different models
AICs <- t(sapply(
	seq(1,length(ffr.fits)),
	function(x,ffr.fits) {
		unlist(lapply(ffr.fits[[x]]$AICs, function(x){mean(unlist(x))}))
	},
	ffr.fits=ffr.fits
))

# name the rows
rownames(AICs) <- labels

# given them slightly more interpretable names
colnames(AICs) <- c(
	"H1",
	"H2.SS",
	"H2.SG",
	"H2.GS",
	"H2.GG",
	"H2.HHI",
	"H2.HHE"
)
AICs <- as.data.frame(AICs)

# break a tie between the two hybrid-hybrid models
AICs[,"H2.HH"] <- pmin(AICs[,"H2.HHI"],AICs[,"H2.HHE"])
AICs[,"H2.HHI"] <- AICs[,"H2.HHE"] <- NULL

# build a table of delta AIC
minAICs <- apply(AICs, 1, min)
dAICs <- AICs - minAICs
#dAICs[dAICs<2] <- 0

# build a table of AIC ranks
rnkAICs <- t(apply(dAICs, 1, rank, ties.method='first'))
colnames(rnkAICs) <- colnames(AICs)
