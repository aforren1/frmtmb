.libPaths(c('C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb.latent))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: hmm
### Title: Hidden Markov models
### Aliases: hmm

### ** Examples

set.seed(11)
n_seq <- 20; len <- 25
G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
mu <- c(0, 3); sg <- c(0.6, 0.6)
dd <- do.call(rbind, lapply(seq_len(n_seq), function(id) {
  s <- integer(len); s[1] <- 1L
  for (t in 2:len) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
  data.frame(id = id, t = seq_len(len),
             y = rnorm(len, mu[s], sg[s]), state = s)
}))
fit <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
fixef(fit)

# smoothed state probabilities and the MAP path
head(hmm_probs(fit))
mean(hmm_viterbi(fit) == dd$state)

# fitted() is the occupancy-weighted mean, not state 1's
cor(fitted(fit), dd$y)

## No test: 
# one state's mean takes its own predictor, random effects included
frm(bf(y ~ 1, mu2 ~ 1 + (1 | id)),
    family = hmm(K = 2, gaussian(), time = t, group = id),
    data = dd)

# covariate-dependent transitions: trans = sets every cell's default
dd$x <- rnorm(nrow(dd))
frm(bf(y ~ 1),
    family = hmm(K = 2, gaussian(), time = t, group = id,
                 init = "estimated", trans = ~x),
    data = dd)
## End(No test)



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
