
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

# # DEBUG
# ffr.fits <- ffr.fits[1:5]

# profile all of the datasets
ffr.cfs <- lapply(
	ffr.fits,
	function(ffr.fit,modeltypes,parameters){
		ffr.cfs <- list()
		# for each model type
		for(modeltype in modeltypes){
			est <- ffr.fit$estimates[[modeltype]]["50%",parameters,"estimate"]

			if(ffr.fit$estimates[[modeltype]]["n",1,1] != 1){
				method <- 'bootstrap'

				# use central interval equivalent to one SD as bounds
				lb <- ffr.fit$estimates[[modeltype]]["16%",parameters,"estimate"]
				ub <- ffr.fit$estimates[[modeltype]]["84%",parameters,"estimate"]

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
				helper <- try(HelpersMG::SEfromHessian(ffr.fit$fits[[modeltype]]@details$hessian, hessian=TRUE))
				if(!inherits(helper, "try-error")){
					std.err <- helper$SE
					ffr.fit$fits[[modeltype]]@details$hessian <- helper$hessian
					ffr.fit$fits[[modeltype]]@vcov <- solve(helper$hessian)
				}else{
					std.err <- summary(ffr.fit$fits[[modeltype]])@coef[parameters,"Std. Error"]
					std.err[is.na(std.err)] <- 0.01
				}

				# try to get profile confidence intervals
				for(parameter in parameters){
					ccff <- try(confint(
						profile(
							ffr.fit$fits[[modeltype]],
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

			# save the model's confidence intervals
			ffr.cfs[[modeltype]] <- cf
		}
		names(ffr.cfs) <- modeltypes
		ffr.cfs <- list(
			study.info = ffr.fit$study.info,
			AICs = unlist(ffr.fit$AICs),
			profiles = ffr.cfs
		)
		return(ffr.cfs)
	},
	modeltypes=modeltypes,
	parameters=parameters
)

# save the mega container which includes all FR fits
save(ffr.cfs,file='../../../../results/R/OnePredTwoPrey_fits_profiled/ffr.fits.prof.HH.Rdata')
