# Review of lane meme: models without me(), lane library against
# rellib-r15, on a wider set than dev/meme-guard.R.
#   Rscript dev/review-meme-guard.R <lane|base> > dev/review-meme-out/guard-<arm>.txt
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(9)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 6)),
                tt = rep(1:6, 10))
d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(n)
d$y2 <- 0.3 * d$x + 0.5 * d$y + rnorm(n)
d$v <- 0.5 * d$x + rnorm(n)
d$v[c(4, 19, 33, 47)] <- NA
d$sy <- 0.3
d$xo <- sample(1:4, n, TRUE)
d$yo <- cut(d$y + rnorm(n), c(-Inf, -0.5, 0.5, 1.5, 2.5, Inf),
            ordered_result = TRUE)
shapes <- list(
  mi = list(f = bf(y ~ mi(v) + x) + bf(v | mi() ~ z) + set_rescor(FALSE),
            fam = NULL,
            pr = set_prior("normal(0, 1)", class = "b", resp = "y")),
  rescor = list(f = mvbf(bf(y ~ x), bf(y2 ~ x)) + set_rescor(TRUE),
                fam = NULL,
                pr = set_prior("lkj(2)", class = "rescor")),
  corsd = list(f = bf(y ~ x + (1 + x | g)), fam = gaussian(),
               pr = set_prior("exponential(1)", class = "sd") +
                 set_prior("lkj(2)", class = "cor") +
                 set_prior("normal(0, 1)", class = "b")),
  ar = list(f = bf(y ~ x + ar(time = tt, gr = g, p = 1, cov = TRUE)),
            fam = gaussian(),
            pr = set_prior("normal(0, 1)", class = "b")),
  ordinal = list(f = bf(yo ~ x), fam = cumulative(threshold = "equidistant"),
                 pr = set_prior("normal(0, 1)", class = "b")),
  mo = list(f = bf(y ~ mo(xo) + x), fam = gaussian(),
            pr = set_prior("normal(0, 1)", class = "b")),
  misdy = list(f = bf(y | mi(sy) ~ x), fam = gaussian(),
               pr = set_prior("normal(0, 1)", class = "b")))
r15 <- function(v) sprintf("%.15g", v)
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("\n==== ", nm, "\n", sep = "")
  dp <- tryCatch(as.data.frame(default_prior(s$f, data = d, family = s$fam)),
                 error = function(e) conditionMessage(e))
  print(dp)
  ds <- tryCatch(as.data.frame(default_prior(s$f, data = d, family = s$fam,
                                             route = "sample")),
                 error = function(e) conditionMessage(e))
  print(ds)
  print(tryCatch(as.data.frame(validate_prior(s$pr, s$f, data = d,
                                              family = s$fam)),
                 error = function(e) conditionMessage(e)))
  fit <- tryCatch(suppressWarnings(frm(s$f, data = d, family = s$fam,
                                       prior = s$pr)),
                  error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    cat("FIT ERROR", fit, "\n")
    next
  }
  p <- fit$obj$env$last.par.best
  cat("par:", r15(p), "\n", fill = 100)
  cat("logLik:", r15(as.numeric(logLik(fit))), "\n")
  p1 <- p + 0.01 * seq_along(p) / length(p)
  ri <- fit$obj$env$random
  cat("f(p1):", r15(fit$obj$fn(if (length(ri)) p1[-ri] else p1)), "\n")
  cat("joint f(p1):", r15(fit$obj$env$f(p1)), "\n")
  print(tryCatch(as.data.frame(prior_summary(fit)),
                 error = function(e) conditionMessage(e)))
  msgs <- character(0)
  dr <- tryCatch(withCallingHandlers(
    suppressWarnings(frm_sample(s$f, data = d, family = s$fam,
                                prior = s$pr, sample_prior = "yes",
                                chains = 1, iter = 300, seed = 4,
                                refresh = 0)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }), error = function(e) conditionMessage(e))
  cat("messages:\n", msgs, "\n", sep = "")
  if (is.character(dr)) {
    cat("SAMPLE ERROR", dr, "\n")
    next
  }
  m <- as.matrix(dr)
  cat("draws cols:", colnames(m), "\n", fill = 100)
  cat("draws digest (1e-8):", digest::digest(round(unname(m), 8)), "\n")
  cat("prior draws:", names(prior_draws(dr)), "\n")
  print(as.data.frame(prior_summary(dr)))
}
cat("DONE\n")
