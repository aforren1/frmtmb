# The one model both changed lanes of the 0.78.0 round reach, meme
# (brms's me() prior classes) and rescorpp (the joint predictive draw
# and per-response density of set_rescor(TRUE) draws): mvbf() with
# rescor and me() terms. Data seed 29 (rm_data() of the test file),
# sampler seed 5, 2 chains of 800 (400 warmup).
#
#   Rscript dev/rel078-overlap.R [lib] > dev/rel078-log/overlap-<n>.txt
#
# lib: a library in front of rellib-r15 and the user library (default
# rellib-r16, the release build). It prints brms 2.23.0's
# default_prior() rows on the formula (no compile), frmtmb's
# sampling-route table, and on the draws: the variables, the prior
# draws' names, the shape of posterior_predict(), the per-draw
# correlation of the predictive residuals against rescor, and
# log_lik(resp =) against a hand density. The test file
# extensions/frmtmb.sample/tests/testthat/test-rescor-me.R
# (dev/rel078-fix-rescorme-test.R) asserts what this prints.
lib <- commandArgs(TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/rellib-r16"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r15",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(NOT_CRAN = "true")
cat("brms", format(packageVersion("brms")), "\n")
src <- "C:/Users/adf44/source/r/frmtmb-wt-rel078/dev/rel078-fix-rescorme-test.R"
ex <- parse(src)
keep <- vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("rm_data", "rm_form", "rm_prior",
                                "rm_quiet")
}, NA)
suppressPackageStartupMessages(library(brms))
for (e in ex[keep]) eval(e)
d <- rm_data()
bf_b <- brms::mvbf(brms::bf(y ~ me(xo, sx)),
                   brms::bf(y2 ~ me(xo, sx) + me(zo, sz))) +
  brms::set_rescor(TRUE)
bp <- brms::default_prior(bf_b, data = d)
cat("== brms 2.23.0 default_prior(), me() and rescor rows\n")
bpm <- as.data.frame(bp)[bp$class %in% c("meanme", "sdme", "corme",
                                          "rescor", "Lme", "Lrescor"),
                         c("prior", "class", "coef", "group", "resp",
                           "lb", "ub")]
print(bpm, row.names = FALSE)
detach("package:brms", unload = FALSE)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    find.package("frmtmb.sample"), format(packageVersion("frmtmb.sample")),
    "\n")
for (e in ex[keep]) eval(e)
q <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  if (is.character(r) && length(r) == 1L) cat("  ", r, "\n") else print(r)
  invisible(r)
}
q("frmtmb default_prior(route = \"sample\"), me() and rescor rows", {
  tab <- default_prior(rm_form(), data = d, route = "sample")
  as.data.frame(tab)[tab$class %in% c("meanme", "sdme", "corme", "rescor"),
                     c("prior", "class", "coef", "group", "resp", "lb",
                       "ub")]
})
ds <- q("frm_sample()", rm_quiet(frm_sample(
  rm_form(), data = d, family = gaussian(), prior = rm_prior(),
  sample_prior = "yes", chains = 2, iter = 800, seed = 5, refresh = 0)))
if (inherits(ds, "frmtmb_draws")) {
  m <- as.matrix(ds)
  q("variables", colnames(m))
  q("prior_draws() names", names(prior_draws(ds)))
  pp <- posterior_predict(ds)
  q("dim(posterior_predict()) and its responses",
    list(dim(pp), dimnames(pp)[[3L]]))
  ep <- posterior_epred(ds)
  if (length(dim(pp)) == 3L) {
    r <- vapply(seq_len(nrow(m)), function(s) {
      stats::cor(pp[s, , 1] - ep[s, , 1], pp[s, , 2] - ep[s, , 2])
    }, 0)
    cat(sprintf("mean per-draw sample correlation %.4f\n", mean(r)))
    # the unbiased reading: the mean product of the standardized
    # residuals, whose expectation is the draw's rescor
    u <- rowMeans((pp[, , 1] - ep[, , 1]) / m[, "sigma_y"] *
                    (pp[, , 2] - ep[, , 2]) / m[, "sigma_y2"])
    dlt <- u - m[, "rescor__y__y2"]
    cat(sprintf(paste("mean standardized product %.4f, mean rescor %.4f,",
                      "z %.3f, z of rescor against 0 %.1f\n"),
                mean(u), mean(m[, "rescor__y__y2"]),
                mean(dlt) / (sd(dlt) / sqrt(length(dlt))),
                mean(m[, "rescor__y__y2"]) / (sd(dlt) / sqrt(length(dlt)))))
  }
  q("log_lik(resp = \"y\") against dnorm(), max relative difference", {
    ll1 <- log_lik(ds, resp = "y")
    ref1 <- dnorm(matrix(d$y, nrow(m), nrow(d), byrow = TRUE),
                  ep[, , 1], m[, "sigma_y"], log = TRUE)
    max(abs(ll1 - ref1)) / max(abs(ref1))
  })
  q("log_lik() against the bivariate normal, max relative difference", {
    rho <- m[, "rescor__y__y2"]
    z1 <- (matrix(d$y, nrow(m), nrow(d), byrow = TRUE) - ep[, , 1]) /
      m[, "sigma_y"]
    z2 <- (matrix(d$y2, nrow(m), nrow(d), byrow = TRUE) - ep[, , 2]) /
      m[, "sigma_y2"]
    refj <- -log(2 * pi) - log(m[, "sigma_y"]) - log(m[, "sigma_y2"]) -
      0.5 * log(1 - rho^2) -
      (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * (1 - rho^2))
    ll <- log_lik(ds)
    c(joint = max(abs(ll - refj)) / max(abs(refj)),
      margins_sum = max(abs(ll - log_lik(ds, resp = "y") -
                              log_lik(ds, resp = "y2"))) / max(abs(refj)))
  })
  q("summary of the me() and rescor variables",
    round(t(apply(m[, grep("^(meanme|sdme|corme|rescor)", colnames(m))], 2,
                  function(v) c(mean = mean(v), sd = sd(v)))), 4))
}
cat("DONE\n")
