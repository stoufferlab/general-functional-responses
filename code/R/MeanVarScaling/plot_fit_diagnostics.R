# Plot the observed vs. fitted for each dataset and model

source('lib/plot_coefs.R')


type <- 'One_Predator_One_Prey'
dir <- paste0('../../../results/R/MeanVarScaling/IndividualFits/', type)
ffr.fits <- bundle_fits(dir)


pdf(paste0('../../../results/R/MeanVarScaling/Residual-vs-Predicted_', type, '.pdf'), 
    height = 8, width = 12)
nModels <- length(ffr.fits[[1]]$fits)
par(mfrow = c(5, nModels), mgp = c(1, 0.2, 0), 
    mar = c(3, 3, 1.5, 0.5), pty = 's',
    cex = 0.7, cex.main = 0.7, tcl = -0.2)
for(i in seq_along(ffr.fits)){
  fits <- ffr.fits[[i]]$fits
  for(m in seq_along(fits)){
    resids <- residuals(fits[[m]])
    pred <- predict(fits[[m]])
    # obs <- fits[[m]]$model[, grep('.var$', names(fits[[m]]$model))]
    wts <- fits[[m]]$model[, grep('weights', names(fits[[m]]$model))]
    plot(pred,
         resids,
         xlab = 'Predicted',
         ylab = 'Residuals', 
         main = paste0(ffr.fits[[i]]$study.info$datasetName, '\n',
                       names(fits)[m]))
    abline(h = 0, lty = 2)
    r <- lm(resids~pred, weights = wts)
    if(length(coef(r)) == r$rank){ # to avoid warning for mean-only model
      lines(range(pred), predict(r, data.frame(pred=range(pred))))
    }
  }
}
dev.off()


###############################
type <- 'One_Predator_Two_Prey'
dir <- paste0('../../../results/R/MeanVarScaling/IndividualFits/', type)
ffr.fits <- bundle_fits(dir)


pdf(paste0('../../../results/R/MeanVarScaling/Residual-vs-Predicted_', type, '.pdf'), 
    height = 8, width = 12)
nModels <- length(ffr.fits[[1]]$fits)
par(mfrow = c(5, nModels), mgp = c(1, 0.2, 0), 
    mar = c(3, 3, 1.5, 0.5), pty = 's',
    cex = 0.7, cex.main = 0.7, tcl = -0.2)
for(i in seq_along(ffr.fits)){
  fits <- ffr.fits[[i]]$fits
  for(p in c(1,2)){
    for(m in seq_along(fits)){
      resids <- residuals(fits[[m]][[p]])
      pred <- predict(fits[[m]][[p]], type = 'response')
      # obs <- fits[[m]][[p]]$model[, grep('.var$', names(fits[[m]][[p]]$model))]
      wts <- fits[[m]][[p]]$model[, grep('weights', names(fits[[m]][[p]]$model))]
      plot(pred,
           resids,
           xlab = 'Predicted', 
           ylab = 'Residuals',
           main = paste0(ffr.fits[[i]]$study.info$datasetName, 
                         ' Prey ', p,
                         '\n',
                         names(fits)[m]))
      abline(h = 0, lty = 2)
      r <- lm(resids~pred, weights = wts)
      if(length(coef(r)) == r$rank){ # to avoid warning for mean-only model
        lines(range(pred), predict(r, data.frame(pred=range(pred))))
      }
    }
  }
}
dev.off()
