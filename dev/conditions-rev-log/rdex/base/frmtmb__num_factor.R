.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: num_factor
### Title: Factor with numeric-coded levels for coordinate covariance
###   structures
### Aliases: num_factor

### ** Examples

# unequally spaced observation times, kept as distances
num_factor(c(0, 1.5, 4))

# planar coordinates for a spatial covariance
levels(num_factor(rep(1:3, 3), rep(1:3, each = 3)))

# ou() reads the distances out of the level labels; a plain factor
# would only give it an ordering
set.seed(1)
tim <- c(0, 1, 1.5, 3)
n_g <- 40
S <- 0.9^2 * exp(-1.2 * abs(outer(tim, tim, "-")))
u <- matrix(rnorm(n_g * length(tim)), n_g) %*% chol(S)
dd <- data.frame(
  y = 1 + as.vector(t(u)) + rnorm(n_g * length(tim), 0, 0.4),
  g = factor(rep(seq_len(n_g), each = length(tim))),
  tim = num_factor(rep(tim, n_g))
)
fit <- frm(bf(y ~ 1 + ou(tim + 0 | g)) + gaussian(), data = dd)
round(VarCorr(fit)$g$cor[, "Estimate", ], 3)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
