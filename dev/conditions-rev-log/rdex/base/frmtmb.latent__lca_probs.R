.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.latent))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: lca_probs
### Title: Posterior class membership of an 'lca()' fit
### Aliases: lca_probs

### ** Examples

set.seed(3)
n <- 200
cl <- rbinom(n, 1, 0.5) + 1
pr <- rbind(c(0.9, 0.85, 0.8, 0.9), c(0.1, 0.15, 0.2, 0.1))
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(row = seq_len(n))
dd$Y <- Y
fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)

p <- lca_probs(fit)
head(p)
attr(p, "entropy")            # classification quality, 0 to 1
table(max.col(p), truth = cl) # the modal assignment, up to relabeling



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
