.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.latent))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: hmm_starts
### Title: Refits of an 'hmm()' model from jittered starting values
### Aliases: hmm_starts

### ** Examples

set.seed(4)
n_seq <- 12; len <- 20
G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
mu <- c(0, 3); sg <- c(0.6, 0.6)
dd <- do.call(rbind, lapply(seq_len(n_seq), function(id) {
  s <- integer(len); s[1] <- 1L
  for (t in 2:len) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
  data.frame(id = id, t = seq_len(len), y = rnorm(len, mu[s], sg[s]))
}))
fit <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
ms <- hmm_starts(fit, n = 3, seed = 1)
ms
logLik(ms$best)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
