.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.latent))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: hmm_viterbi
### Title: Most likely state path of an hmm fit (Viterbi)
### Aliases: hmm_viterbi

### ** Examples

set.seed(12)
dd <- data.frame(id = 1, t = 1:120)
s <- integer(120); s[1] <- 1L
G <- matrix(c(0.92, 0.08, 0.15, 0.85), 2, 2, byrow = TRUE)
for (t in 2:120) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
dd$y <- rnorm(120, c(0, 3)[s], 0.5)
fit <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
v <- hmm_viterbi(fit)
table(v, truth = s)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
