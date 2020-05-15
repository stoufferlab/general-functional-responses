
# load necessary libraries
library(bbmle)
library(HelpersMG)

# make sure we have functions locally
source('../../lib/holling_method_one_predator_two_prey.R')

# only profile the hybrid model
modeltypes <- c(
	"Holling II Hybrid Hybrid identity",
	"Holling II Hybrid Hybrid exp"
)

# only get confidence intervals on the phis
# everything else takes way too long
parameters <- c("phi_ij","phi_ji")

# read in unprofiled fits
ffr.fits <- readRDS(
	file='../../../../results/R/OnePredTwoPrey_ffr.fits.Rdata'
)

# # scrape out the AIC values for the different models
# AICs <- t(sapply(
# 	seq(1,length(ffr.fits)),
# 	function(x,ffr.fits) {
# 		unlist(lapply(ffr.fits[[x]]$fits, AIC))
# 	},
# 	ffr.fits=ffr.fits
# ))
# colnames(AICs) <- c(
# 	"H1",
# 	"H2.SS",
# 	"H2.SG",
# 	"H2.GS",
# 	"H2.GG",
# 	"H2.HHI",
# 	"H2.HHE"
# )
# AICs <- as.data.frame(AICs)

# # figure out which version of the hybrid-hybrid model fit best
# better.model <- apply(
# 	AICs,
# 	MAR=1,
# 	function(x){
# 		ifelse(
# 			x["H2.HHI"]<x["H2.HHE"],
# 			modeltypes[1],
# 			modeltypes[2]
# 		)
# 	}
# )

# container for profiled CIs of all datasets
ffr.cfs <- list()

# profile each dataset
for(x in 1:length(ffr.fits)){
	ffr.fit <- ffr.fits[[x]]
	cfs <- list()
	for(model in modeltypes){
	# model <- better.model[x]

	est <- ffr.fit$estimates[[model]]["50%",parameters,"estimate"]

	if(ffr.fit$estimates[[model]]["n",1,1] != 1){
		method <- 'bootstrap'

		# use central interval equivalent to one SD as bounds
		lb <- ffr.fit$estimates[[model]]["16%",parameters,"estimate"]
		ub <- ffr.fit$estimates[[model]]["84%",parameters,"estimate"]

		# create a data frame with the confidence intervals
		cf <- data.frame(
			row.names=parameters,
			method=rep("bootstrap",length(parameters)),
			lb=lb,
			est=est,
			ub=ub,
			stringsAsFactors = FALSE
		)
	}else{
		cf <- data.frame(
			row.names=parameters,
			method=rep("foobar",length(parameters)),
			lb=rep(NA,length(parameters)),
			est=est,
			ub=rep(NA,length(parameters)),
			stringsAsFactors = FALSE
		)

		# make sure we have a well-behaved Hessian to use to jumpstart profiling
		helper <- try(HelpersMG::SEfromHessian(ffr.fit$fits[[model]]@details$hessian, hessian=TRUE))
		if(!inherits(helper, "try-error")){
			std.err <- helper$SE
			ffr.fit$fits[[model]]@details$hessian <- helper$hessian
			ffr.fit$fits[[model]]@vcov <- solve(helper$hessian)
		}else{
			std.err <- summary(ffr.fit$fits[[model]])@coef[parameters,"Std. Error"]
			std.err[is.na(std.err)] <- 0.01
		}

		# try to get profile confidence intervals
		for(parameter in parameters){
			ccff <- try(confint(
				profile(
					ffr.fit$fits[[model]],
					which=parameter,
					try_harder=TRUE,
					tol.newmin=Inf,
					std.err=std.err[parameter],
					trace=TRUE,
					zmax=2.5
				),
				level=0.68
			))
			if(!inherits(cf, "try-error")){
				cf[parameter,"method"] <- "profile"
				cf[parameter,"lb"] <- ccff[1]
				cf[parameter,"ub"] <- ccff[2]
			}
		}
	}

	if(grepl("exp",model)){
		cf[,c("lb","est","ub")] <- exp(cf[,c("lb","est","ub")])
	}

	cfs[[model]] <- cf
	}

	ffr.cfs[[x]] <- list(
		study.info = ffr.fit$study.info,
		AICs = unlist(lapply(ffr.fit$AICs, function(x){mean(unlist(x))})),
		profiles = cfs
	)
}

# save the mega container which includes confidence intervals of all datasets
save(ffr.cfs,file='../../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata')
