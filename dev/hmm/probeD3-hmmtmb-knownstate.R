# Probe D3: the hmmTMB sharp edge found while building probe D2.
#
# hmmTMB reads a data column named `state` as KNOWN STATES and silently
# conditions the likelihood on it. A simulation study that keeps the
# generating state sequence in the data frame - the obvious thing to do -
# therefore compares against a DIFFERENT likelihood, with no message.
# The three fits below are identical except for which columns the data
# frame carries.
#
# Run: Rscript dev/hmm/probeD3-hmmtmb-knownstate.R

lib <- Sys.getenv("FRMTMB_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(hmmTMB))
source("dev/hmm/hmm-common.R")

K <- 2L
set.seed(2026)
N <- 25L
Tg <- 30L
n <- N * Tg
G_true <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
mu_true <- c(0, 3)
sigma_true <- c(0.6, 0.6)
sd_b <- c(0.7, 0.5)
b_true <- cbind(rnorm(N, 0, sd_b[1]), rnorm(N, 0, sd_b[2]))
dat <- do.call(rbind, lapply(seq_len(N), function(g) {
  s <- integer(Tg)
  s[1] <- sample.int(K, 1, prob = stat_dist(G_true))
  for (t in seq_len(Tg - 1L)) {
    s[t + 1L] <- sample.int(K, 1, prob = G_true[s[t], ])
  }
  data.frame(ID = g, t = seq_len(Tg), state = s,
             y = rnorm(Tg, mu_true[s] + b_true[g, s], sigma_true[s]))
}))

run <- function(d, tag) {
  hid <- MarkovChain$new(data = d, n_states = 2,
                         initial_state = "stationary")
  obs <- Observation$new(data = d, dists = list(y = "norm"),
                         n_states = 2,
                         par = list(y = list(mean = mu_true,
                                             sd = sigma_true)),
                         formulas = list(y = list(mean = ~ 1, sd = ~ 1)))
  hm <- HMM$new(obs = obs, hid = hid)
  suppressWarnings(hm$fit(silent = TRUE))
  ph <- hm$obs()$par()[, , 1]
  Gh <- hm$hid()$tpm()[, , 1]
  rbg <- split(seq_len(nrow(d)), d$ID)
  # the ordinary (states unknown) forward algorithm at hmmTMB's own
  # fitted parameters
  ours <- sum(vapply(rbg, function(r)
    fwd_num_tv(lpmat_gauss(d$y, ph["y.mean", ], ph["y.sd", ]),
               function(rr) Gh, r, stat_dist(Gh)), numeric(1)))
  cat(sprintf("  %-28s hmmTMB llk %14.6f   plain forward at the same\n",
              tag, hm$llk()))
  cat(sprintf("  %-28s parameters %14.6f   difference %+11.6f\n", "",
              ours, ours - hm$llk()))
  invisible(hm)
}

cat("== probe D3: hmmTMB and a column named `state` ==\n\n")
run(dat, "columns ID, t, state, y")
run(dat[, c("ID", "t", "y")], "columns ID, t, y")
run(dat[, c("ID", "y")], "columns ID, y")
cat("\nThe first fit maximizes the COMPLETE-DATA likelihood (states\n",
    "known); the other two maximize the marginal HMM likelihood and\n",
    "reproduce the plain forward algorithm to the last digit.\n")
