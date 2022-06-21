# Get the p-value for whether the specified coefficient of a fitted model
# differs from a specified value (e.g., val = 1) rather than zero.
# (From: https://stats.stackexchange.com/questions/111559/test-model-coefficient-regression-slope-against-some-value)
ttest <- function(fit, coefnum, val = 1){
  co <- coef(summary(fit))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  2 * pt(abs(tstat), fit$df.residual, lower.tail = FALSE)
}