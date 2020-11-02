
# we need the bundle_fits function
source('../../lib/plot_coefs.R')

# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredTwoPrey_fits')

# scrape out the dataset names
labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))

# scrape out the AIC values for the different models and convert them to BICs
BICs <- t(sapply(
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
		# go through the AIC values and correct to make them BIC values
		BIC <- unlist(sapply(
			1:length(k),
			function(x,AICs,k,SS){
				mean(unlist(AICs[[x]])+(log(SS)-2)*k[[x]])
			},
			AICs=ffr.fit$AICs,
			k=k,
			SS=SS
		))
		return(BIC)
	},
	ffr.fits=ffr.fits
))

# name the rows
rownames(BICs) <- labels

# given them slightly more interpretable names
colnames(BICs) <- c(
	"H1",
	"H2.SS",
	"H2.SG",
	"H2.GS",
	"H2.GG",
	"H2.HHI",
	"H2.HHE"
)
BICs <- as.data.frame(BICs)

# break a tie between the two hybrid-hybrid models
BICs[,"H2.HH"] <- pmin(BICs[,"H2.HHI"],BICs[,"H2.HHE"])
BICs[,"H2.HHI"] <- BICs[,"H2.HHE"] <- NULL

# build a table of delta BIC
minBICs <- apply(BICs, 1, min)
dBICs <- BICs - minBICs

# build a table of BIC ranks
rnkBICs <- t(apply(dBICs, 1, rank, ties.method='first'))
colnames(rnkBICs) <- colnames(BICs)

# print out the BIC ranked info
# remove the SG and GS models
dBICs <- dBICs[,c("H1","H2.SS","H2.GG","H2.HH")]
rnkBICs <- rnkBICs[,c("H1","H2.SS","H2.GG","H2.HH")]

# recompute the delta BICs
dBICs <- dBICs - apply(dBICs,1,min)

# rerank this subset of models
rnkBICs <- t(apply(rnkBICs,1,rank))

# Define delta BICc cut-off for "indistinguishably well performing" models
delBICcutoff <- 2
