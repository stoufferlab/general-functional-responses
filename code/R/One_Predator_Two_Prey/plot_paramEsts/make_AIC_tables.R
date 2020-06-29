
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

# print out the AIC ranked info
# remove the SG and GS models
dAICs <- dAICs[,c("H1","H2.SS","H2.GG","H2.HH")]
rnkAICs <- rnkAICs[,c("H1","H2.SS","H2.GG","H2.HH")]

# recompute the delta AICs
dAICs <- dAICs - apply(dAICs,1,min)

# rerank this subset of models
rnkAICs <- t(apply(rnkAICs,1,rank))

# Define delta AICc cut-off for "indistinguishably well performing" models
delAICcutoff <- 2

# statistics about when the phi model ranks first or ties for first
message("model ranks")
message(paste(
  "phi itself:",
  sum(rnkAICs[,"H2.HH"]==1),
  round(sum(rnkAICs[,"H2.HH"]==1)/nrow(rnkAICs)*100),
  "%"
))
message(paste(
  "phi tied:",
  sum(dAICs[,"H2.HH"]<=2)-sum(rnkAICs[,"H2.HH"]==1),
  round((sum(dAICs[,"H2.HH"]<=2)-sum(rnkAICs[,"H2.HH"]==1))/nrow(rnkAICs)*100),
  "%"
))