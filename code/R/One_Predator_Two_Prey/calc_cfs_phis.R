
modeltype <- "Holling II Hybrid Hybrid"
parameters <- c("phi_ij","phi_ji")

ffr.cfs <- lapply(
	ffr.fits,
	function(x,modeltype,parameters){
		if(x$estimates[[modeltype]]["n",1,1] != 1){
			return(NULL)
		}else{
			foobar <- list()
			foobar[[modeltype]] <- try(confint(
				x$fits[[modeltype]],
				parm=parameters,
				try_harder=TRUE,
				level=0.68,
				tol.newmin=Inf,
				quietly=FALSE
			))
			foobar
		}
	},
	modeltype=modeltype,
	parameters=parameters
)

# save the mega container which includes all FR fits
save(ffr.cfs,file='../../../results/R/ffr.cfs_OnePredTwoPrey.Rdata')
