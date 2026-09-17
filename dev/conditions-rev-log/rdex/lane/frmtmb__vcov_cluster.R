.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: vcov_cluster
### Title: Cluster-robust (sandwich) covariance
### Aliases: vcov_cluster

### ** Examples

set.seed(1)
G <- 40
dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                 x = rnorm(G * 6))
# cluster-specific error scale the model does not describe
dd$y <- 1 + 0.5 * dd$x +
  rnorm(G * 6, 0, rep(runif(G, 0.3, 2.5), each = 6))
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           REML = FALSE)

sqrt(diag(vcov(fit)))                          # model-based
sqrt(diag(vcov_cluster(fit, ~ g, "CR1")))      # cluster-robust

# feeds straight into the inference methods, which need the whole
# outer parameter vector and pick up the t(G - 1) reference from it
confint(fit, parm = "x",
        vcov = vcov_cluster(fit, ~ g, "CR1", full = TRUE))



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
