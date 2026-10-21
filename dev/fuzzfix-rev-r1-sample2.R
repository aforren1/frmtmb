# Reviewer, fix round 1: frm_sample() formula route on models whose
# posterior is proper although sum(w) > n (data of dev/fuzzfix-rev-wu.R).
#   Rscript dev/fuzzfix-rev-r1-sample2.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), w = runif(n, 0.5, 2))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, 0.8)
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
run <- function(lbl, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-40s ERROR %s\n", lbl, substr(conditionMessage(r), 1, 120)))
  } else {
    dm <- as.matrix(r)
    cn <- grep("^shape|^sigma|betad", colnames(dm), value = TRUE)[1]
    rh <- max(posterior::summarise_draws(r)$rhat, na.rm = TRUE)
    cat(sprintf("%-40s ok %s median %.4g max rhat %.3f\n", lbl, cn,
                stats::median(dm[, cn]), rh))
  }
}
run("Gamma, default priors (gamma on shape)",
    frm_sample(bf(yg | weights(w) ~ x + (1 | obs)), data = d,
               family = Gamma("log"), chains = 2, iter = 1000, seed = 1,
               refresh = 0))
run("gaussian, lognormal(0, 1) on sigma",
    frm_sample(bf(y | weights(w) ~ x + (1 | obs)), data = d,
               family = gaussian(), chains = 2, iter = 1000, seed = 1,
               refresh = 0, prior = set_prior("lognormal(0, 1)",
                                              class = "sigma")))
