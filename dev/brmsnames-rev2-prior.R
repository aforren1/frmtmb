## Reviewer recheck: why does frm_sample(prior = set_prior(class = "b"))
## give different draws between arms when the flat run is identical?
##   Rscript dev/brmsnames-rev2-prior.R base|lane
## Data seed 41 (the gauss model of dev/brmsnames-rev2-natural.R).
arm <- commandArgs(trailingOnly = TRUE)[1L]
libs <- list(
  base = c("C:/Users/adf44/source/r/rellib-r3"),
  lane = c("C:/Users/adf44/source/r/brmsnames-lib",
           "C:/Users/adf44/source/r/rellib-r3"))
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
set.seed(41)
n <- 160; G <- 8
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), length.out = n)))
u <- rnorm(G, 0, 0.7)
d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(n, 0, 4)
fit <- q(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d))
pri <- set_prior("normal(0, 1)", class = "b")
ns <- asNamespace("frmtmb.sample")
defaults <- ns$default_priors_for(fit)
cat("== default priors ==\n"); print(defaults)
res <- ns$sample_resolve_priors(fit, pri, base = defaults)
cat("== resolved entries ==\n")
str(lapply(res$entries %||% res, function(e) e[c("comp", "idx", "scale")]),
    max.level = 2)
msgs <- character(0)
ds <- withCallingHandlers(
  suppressWarnings(frm_sample(fit, prior = pri, chains = 1, iter = 60,
                              refresh = 0, seed = 9)),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
  })
cat("== messages ==\n"); cat(msgs, sep = "")
cat("== prior_summary(ds) ==\n"); print(prior_summary(ds))
# the log density at a fixed parameter vector, from the stanfit's own
# log_prob: the flat part is identical between arms, so a difference is
# the prior
sf <- ds$stanfit
full <- fit$obj$env$last.par.best
cat("npar stanfit", rstan::get_num_upars(sf), " full par", length(full), "\n")
cat("== log_prob at the ML mode and at shifted points ==\n")
for (sh in c(0, 0.3, -0.5)) {
  v <- full
  v[1:4] <- v[1:4] + sh
  cat(sprintf("shift %4.1f: %.10f\n", sh,
              rstan::log_prob(sf, v, adjust_transform = FALSE)))
}
for (j in 1:4) {
  v <- full
  v[j] <- v[j] + 0.5
  cat(sprintf("par %d + 0.5: %.10f\n", j,
              rstan::log_prob(sf, v, adjust_transform = FALSE)))
}
cat("== first 3 draws ==\n")
print(ds$draws[1:3, 1:4])
