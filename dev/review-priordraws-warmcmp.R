# Reviewer: the warmup the fix round rebuilds on demand against the copy
# the previous build stored, on models the lane did not use.
#   Rscript dev/review-priordraws-warmcmp.R <old|new>
# old: the previous frmtmb.sample (wt-priordraws-rv-lib/floor) on the
# lane core; new: the lane library. Saves warmup-<arm>.rds.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "old") "C:/Users/adf44/source/r/wt-priordraws-rv-lib/floor",
            "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1)
d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(6, 0, 0.5)[d$g] + rnorm(60)
d$t <- rep(1:12, 5)
d$ye <- d$y + rexp(60)
pr <- set_prior("normal(0, 1)", class = "b")
run <- function(f, fam = gaussian(), prior = pr, fitfirst = FALSE, ...) {
  if (fitfirst) return(q(frm_sample(q(frm(bf(f), data = d, family = fam)), prior = prior,
    chains = 2, iter = 300, warmup = 200, seed = 13, refresh = 0, ...)))
  q(frm_sample(bf(f), data = d, family = fam, prior = prior, chains = 2,
               iter = 300, warmup = 200, seed = 13, refresh = 0, ...))
}
cases <- list(
  re_cor = function() run(y ~ x + (1 + x | g), sample_prior = "yes"),
  mo_student = function() run(y ~ mo(m) + x, student(), sample_prior = "yes"),
  ar_cov = function() run(y ~ x + ar(time = t, gr = h, cov = TRUE),
                          prior = pr + set_prior("normal(0, 0.5)", class = "ar")),
  gev = function() run(ye ~ x, gen_extreme_value(), fitfirst = TRUE),
  thin = function() run(y ~ x + (1 | g), thin = 2)
)
out <- list()
for (nm in names(cases)) {
  ds <- tryCatch(cases[[nm]](), error = function(e) e)
  if (inherits(ds, "error")) { cat(nm, "ERROR", conditionMessage(ds), "\n"); next }
  a <- tryCatch(unclass(as_draws_array(ds, inc_warmup = TRUE)), error = function(e) conditionMessage(e))
  mc <- tryCatch(dim(as.mcmc(ds, inc_warmup = TRUE)[[1]]), error = function(e) conditionMessage(e))
  rh <- tryCatch(rhat(ds, inc_warmup = TRUE), error = function(e) conditionMessage(e))
  nr <- tryCatch(nrow(as.data.frame(ds, inc_warmup = TRUE)), error = function(e) conditionMessage(e))
  cat(sprintf("%s: stored warmup %s, size %.0f KB, as_draws_array %s, as.mcmc %s, as.data.frame rows %s\n",
              nm, !is.null(ds$warmup), object.size(ds) / 1024,
              paste(dim(a), collapse = "x"), paste(mc, collapse = "x"), nr))
  out[[nm]] <- list(draws = ds$draws, a = a, rh = rh, mc = mc, size = object.size(ds))
}
# several imputations
dm <- d
dm$x[c(3, 17, 40)] <- NA
dl <- list(d, transform(d, y = y + rnorm(60, 0, 0.1)))
mf <- tryCatch(q(frm_multiple(bf(y ~ x + (1 | g)), data = dl, family = gaussian())),
               error = function(e) e)
if (!inherits(mf, "error")) {
  dsm <- tryCatch(q(frm_sample(mf, prior = pr, chains = 2, iter = 300, warmup = 200,
                               seed = 13, refresh = 0)), error = function(e) e)
  if (inherits(dsm, "error")) cat("multiple ERROR", conditionMessage(dsm), "\n") else {
    a <- unclass(as_draws_array(dsm, inc_warmup = TRUE))
    cat("multiple: stored warmup", !is.null(dsm$warmup), " as_draws_array",
        paste(dim(a), collapse = "x"), "\n")
    out$multiple <- list(draws = dsm$draws, a = a)
  }
} else cat("frm_multiple ERROR", conditionMessage(mf), "\n")
saveRDS(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log/warmup-%s.rds", arm))
