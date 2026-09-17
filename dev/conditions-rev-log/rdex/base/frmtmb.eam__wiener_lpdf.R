.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.eam))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: wiener_lpdf
### Title: The Wiener first-passage log density, for another package's
###   family
### Aliases: wiener_lpdf

### ** Examples

# the density over the upper boundary integrates to the probability
# of ever reaching it, which the diffusion has in closed form
v <- 1.2; a <- 1.5; w <- 0.3
t <- seq(1e-5, 12, length.out = 200000)
p_up <- sum(exp(wiener_lpdf(t, v, a, w, 1))) * diff(t)[1]
p_lo <- (exp(-2 * v * a) - exp(-2 * v * a * w)) / (exp(-2 * v * a) - 1)
c(quadrature = p_up, analytic = 1 - p_lo)

# and it differentiates on an RTMB tape, which is the point
tp <- RTMB::MakeTape(function(v) sum(wiener_lpdf(c(0.4, 0.9), v,
                                                 1.5, 0.3, 1)), 1)
tp$jacobian(1.2)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
