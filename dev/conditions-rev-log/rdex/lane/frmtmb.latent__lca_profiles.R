.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.latent))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: lca_profiles
### Title: The fitted item profiles of an 'lca()' fit
### Aliases: lca_profiles

### ** Examples

set.seed(2)
n <- 200
cl <- rbinom(n, 1, 0.5) + 1
pr <- rbind(c(0.9, 0.85, 0.8), c(0.1, 0.15, 0.2))
Y <- matrix(0L, n, 3)
for (j in 1:3) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(row = seq_len(n))
dd$Y <- Y
fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)

pf <- lca_profiles(fit)
pf
pf[[1]]                       # item 1's K x C table
attr(pf, "class_sizes")



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
