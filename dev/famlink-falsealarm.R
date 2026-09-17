# False-alarm measurement for the mean-link refusal (item 3).
#
# Fits EVERY (family, link) pair brms 2.23.0 accepts, on data that family
# is built for, in ONE arm:
#
#   Rscript dev/famlink-falsealarm.R base
#   Rscript dev/famlink-falsealarm.R lane
#
# and writes dev/famlink-falsealarm-<arm>.tsv. dev/famlink-falsealarm-sum.R
# joins the two. A false alarm is a pair the base FITTED and the lane
# refuses; a regression is a pair whose logLik moved between the arms.
# The pair list is read from the lane's generated table, so both arms fit
# the same list. Seed 20260916 for every data set, one data set per
# family, n = 120.
arm <- commandArgs(trailingOnly = TRUE)[1]
lib <- switch(arm,
              base = "C:/Users/adf44/source/r/rellib-r3",
              lane = "C:/Users/adf44/source/r/famlink-lib",
              stop("arm must be base or lane"))
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

# the table is read from the SOURCE file so the base arm, whose frmtmb has
# no table, fits the same pairs
tab_env <- new.env()
sys.source("R/links-brms.R", envir = tab_env)
sets <- tab_env$brms_mu_links

n <- 120
make_data <- function(fam) {
  set.seed(20260916)
  x <- rnorm(n)
  y <- switch(fam,
    gaussian = , student = , skew_normal = , huber = , asym_laplace = ,
    exgaussian = 0.5 + rnorm(n, 0, 0.08),
    lognormal = , shifted_lognormal = rlnorm(n, 1, 0.3),
    hurdle_lognormal = ifelse(runif(n) < 0.2, 0, rlnorm(n, 1, 0.3)),
    zero_inflated_asym_laplace = ifelse(runif(n) < 0.2, 0,
                                        0.5 + rnorm(n, 0, 0.08)),
    Gamma = , weibull = , exponential = , tweedie = , inverse.gaussian = ,
    cox = rgamma(n, 3, 1),
    hurdle_gamma = ifelse(runif(n) < 0.2, 0, rgamma(n, 3, 1)),
    poisson = , negbinomial = , nbinom1 = , geometric = , compois = ,
    zero_inflated_poisson = , zero_inflated_negbinomial = ,
    hurdle_poisson = rpois(n, 3),
    binomial = , beta_binomial = , zero_inflated_binomial = rbinom(n, 10, 0.4),
    bernoulli = rbinom(n, 1, 0.4),
    beta = rbeta(n, 4, 6),
    zero_inflated_beta = ifelse(runif(n) < 0.2, 0, rbeta(n, 4, 6)),
    von_mises = atan2(sin(0.5 + rnorm(n, 0, 0.5)), cos(0.5 + rnorm(n, 0, 0.5))),
    cumulative = , sratio = , cratio = , acat = factor(sample(1:4, n, TRUE),
                                                      ordered = TRUE),
    categorical = factor(sample(c("a", "b", "c"), n, TRUE)),
    stop("no data for ", fam))
  data.frame(y = y, x = x, tr = 10)
}
form_for <- function(fam) {
  if (fam %in% c("binomial", "beta_binomial", "zero_inflated_binomial")) {
    y | trials(tr) ~ x
  } else {
    y ~ x
  }
}
ctor_for <- function(fam) {
  switch(fam, beta = "Beta", fam)
}

rows <- list()
for (fam in names(sets)) {
  if (fam == "multinomial") next   # no link argument to vary
  d <- make_data(fam)
  for (lk in sets[[fam]]) {
    t0 <- proc.time()[["elapsed"]]
    warn <- character(0)
    res <- withCallingHandlers(tryCatch({
      f <- frmtmb:::family_registry[[ctor_for(fam)]](link = lk)
      fit <- suppressMessages(frm(form_for(fam), d, family = f))
      list(status = "fit", ll = as.numeric(logLik(fit)),
           conv = fit$opt$convergence)
    }, error = function(e) {
      list(status = "error", ll = NA_real_, conv = NA_integer_,
           msg = conditionMessage(e))
    }), warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    rows[[length(rows) + 1L]] <- data.frame(
      family = fam, link = lk, status = res$status,
      loglik = sprintf("%.10f", res$ll), convergence = res$conv,
      warnings = length(warn),
      message = substr(gsub("[\r\n\t]+", " ", res$msg %||% ""), 1, 160),
      secs = round(proc.time()[["elapsed"]] - t0, 2))
    cat(sprintf("%-28s %-14s %-6s %s\n", fam, lk, res$status,
                if (identical(res$status, "fit")) sprintf("%.6f", res$ll)
                else substr(res$msg, 1, 80)))
  }
}
out <- do.call(rbind, rows)
write.table(out, paste0("dev/famlink-falsealarm-", arm, ".tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("arm", arm, ":", nrow(out), "pairs,", sum(out$status == "fit"), "fitted\n")
