
library(bbmle)

# we need the bundle_fits function
source('../../lib/plot_coefs.R')

# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredTwoPrey_fits')

# scrape out the dataset names
labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))

# scrape out the AIC values for the different models and convert them to AICc
AICcs <- t(sapply(
	seq(1,length(ffr.fits)),
	function(i,ffr.fits){
		ffr.fit <- ffr.fits[[i]]
		# number of parameters
		k <- lapply(
			ffr.fit$fits,
			function(x){
				length(coef(x))
			}
		)
		# number of observations (times 2 because two sets of observations per sample)
		if("data.Nconsumed1.mean" %in% names(ffr.fit$study.info)){
			SS <- sum(ffr.fit$study.info$data.n)
		}else{
			SS <- length(ffr.fit$study.info$data.Nconsumed1)
		}
		SS <- 2 * SS
		# go through the AIC values and correct to make them AICc values
		AICc <- unlist(sapply(
			1:length(k),
			function(x,AICs,k,SS){
				mean(unlist(AICs[[x]])+(2*(k[[x]])**2 + 2 * k[[x]])/(SS-k[[x]]-1))
			},
			AICs=ffr.fit$AICs,
			k=k,
			SS=SS
		))
		return(AICc)
	},
	ffr.fits=ffr.fits
))

# name the rows
rownames(AICcs) <- labels

# given them slightly more interpretable names
colnames(AICcs) <- c(
	"H1",
	"H2.SS",
	"H2.SG",
	"H2.GS",
	"H2.GG",
	"H2.HHI",
	"H2.HHE"
)
AICcs <- as.data.frame(AICcs)

# break a tie between the two hybrid-hybrid models
AICcs[,"H2.HH"] <- pmin(AICcs[,"H2.HHI"],AICcs[,"H2.HHE"])
AICcs[,"H2.HHI"] <- AICcs[,"H2.HHE"] <- NULL

# build a table of delta AICc
minAICcs <- apply(AICcs, 1, min)
dAICcs <- AICcs - minAICcs

# build a table of AICc ranks
rnkAICcs <- t(apply(dAICcs, 1, rank, ties.method='first'))
colnames(rnkAICcs) <- colnames(AICcs)

# print out the AICc ranked info
# remove the SG and GS models
dAICcs <- dAICcs[,c("H1","H2.SS","H2.GG","H2.HH")]
rnkAICcs <- rnkAICcs[,c("H1","H2.SS","H2.GG","H2.HH")]

# recompute the delta AICcs
dAICcs <- dAICcs - apply(dAICcs,1,min)

# rerank this subset of models
rnkAICcs <- t(apply(rnkAICcs,1,rank))

# Define delta AICcc cut-off for "indistinguishably well performing" models
delAICccutoff <- 2
