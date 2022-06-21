# Plot the estimates from the individually-fitted datasets

source('lib/plot_coefs.R')

empty.space <- function(x){ plot(1,1, type = 'n', ann = FALSE, axes = FALSE) }

pch.factor <- 'Replacement'
bg.color.factor <- 'Clarity'
p.level <- 0.1

type <- 'One_Predator_One_Prey'
dir <- paste0('../../../results/R/MeanVarScaling/IndividualFits/', type)
ffr.fits <- bundle_fits(dir)


model.type = 'lm'

pdf(paste0('../../../results/R/MeanVarScaling/Coefs_', type, '-', model.type,'.pdf'), 
    height = 10, width = 12)

par(mfrow = c(3, 3))
par(mar = c(3, 6, 0.5, 0.5))
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- model.type
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')


empty.space()

#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- paste0(model.type,'.q')
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)


plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- paste0(model.type,'.main')
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)


plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'c' # same order
#~~~~~~~~~~~~~~
# myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Predator effect ', (italic(c)))), line = 1.5)
legend('topright', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomleft', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'c'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Predator effect ', (italic(c)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- paste0(model.type,'.int')
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'd' # same order
#~~~~~~~~~~~~~~
# myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Predator effect ', (italic(d)))), line = 1.5)
legend('topright', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomleft', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'd'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Predator effect ', (italic(d)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')
 
dev.off() 


#########################################

# Note need for specifying xlims

model.type = 'glm'

pdf(paste0('../../../results/R/MeanVarScaling/Coefs_', type, '-', model.type,'.pdf'), 
    height = 10, width = 12)

par(mfrow = c(3, 3))
par(mar = c(3, 6, 0.5, 0.5))
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- model.type
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level)
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

empty.space()
empty.space()

#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- paste0(model.type,'.main')
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level,
           xlim = c(-2,6))
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'c' # same order
#~~~~~~~~~~~~~~
# myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level,
           xlim = c(-5,5))
title(xlab = expression(paste('Predator effect ', (italic(c)))), line = 1.5)
legend('bottomleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'c'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level,
           xlim = c(-5,5))
title(xlab = expression(paste('Predator effect ', (italic(c)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~
model <- paste0(model.type,'.int')
#~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~
parameter <- 'b'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 1, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level,
           xlim = c(-2,6))
title(xlab = expression(paste('Scaling exponent ', (italic(b)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'd' # same order
#~~~~~~~~~~~~~~
# myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level,
           xlim = c(-0.5,0.5))
title(xlab = expression(paste('Predator effect ', (italic(d)))), line = 1.5)
legend('bottomleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

#~~~~~~~~~~~~~~
parameter <- 'd'
#~~~~~~~~~~~~~~
myorder <- order.of.fits(ffr.fits, model, parameter, type)

plot.coefs(ffr.fits[myorder], type, model, parameter, vertLines = 0, 
           bg.color.factor = bg.color.factor, pch.factor = pch.factor, 
           p.level = p.level,
           xlim = c(-0.5,0.5))
title(xlab = expression(paste('Predator effect ', (italic(d)))), line = 1.5)
legend('topleft', 
       title = pch.factor,
       legend = c("True", "False"), 
       pch = c(23, 21),
       bty = 'n')
legend('bottomright', 
       legend = bquote(italic(p) < .(p.level)), 
       pch = 22,
       pt.bg = c('black'),
       bty = 'n')

dev.off() 

  
  