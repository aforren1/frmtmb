# Reviewer of lane ordmix: false alarms of the flat-direction warnings,
# and the fragility of ordinal-mixture fits, on identified models.
# R replicates per model, seeds 1..R, n = 500. Every warning is
# recorded (none suppressed), with errors, the convergence code and
# whether the gradient at the optimum is finite.
# Usage: Rscript dev/ordmix-rev-falsealarm.R <model> [R]
# Output: dev/ordmix-rev-fa-log/<model>.txt (driver redirects stdout)
args <- commandArgs(TRUE)
model <- args[1]
R <- if (length(args) >= 2) as.integer(args[2]) else 20L
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("model", model, "R", R, "lib", find.package("frmtmb"), "\n")
cut4 <- function(lat) 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
gen <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat_none <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  lat_mu <- ifelse(cls == 1, 2 * x, -1.5 * x) + rlogis(n)
  lat_pr <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rnorm(n)
  y <- cut4(lat_none)
  # hurdle: class 1 has hu 0.1, class 2 hu 0.4, and hu moves with z
  hu <- plogis(ifelse(cls == 1, -2, -0.4) + 0.8 * z)
  yh <- ifelse(runif(n) < hu, 0L, y)
  # order = "mu", no mu predictor, theta moving with z, a cumulative
  # and an sratio component at one threshold vector (-1, 0, 1)
  tau <- c(-1, 0, 1)
  pc <- diff(c(0, plogis(tau), 1))
  h <- plogis(tau)
  ps <- c(h[1], (1 - h[1]) * h[2], (1 - h[1]) * (1 - h[2]) * h[3],
          (1 - h[1]) * (1 - h[2]) * (1 - h[3]))
  th <- plogis(0.3 + 1.5 * z)
  c1 <- runif(n) < th
  ysh <- vapply(seq_len(n), function(i) {
    sample.int(4L, 1L, prob = if (c1[i]) pc else ps)
  }, 1L)
  data.frame(x, z, cls, y, yh, ysh, ymu = cut4(lat_mu), ypr = cut4(lat_pr),
             gr = factor(sample(letters[1:2], n, TRUE)))
}
spec <- switch(model,
  cum2 = list(f = bf(y ~ x), fam = mixture(cumulative(), cumulative())),
  probit2 = list(f = bf(ypr ~ x),
                 fam = mixture(cumulative("probit"), cumulative("probit"))),
  cum_sr_clog = list(f = bf(y ~ x),
                     fam = mixture(cumulative(), sratio("cloglog"))),
  mu2 = list(f = bf(ymu ~ x),
             fam = mixture(cumulative(), cumulative(), order = "mu")),
  hurdle_huz = list(f = bf(yh ~ x, hu1 ~ z, hu2 ~ z),
                    fam = mixture(hurdle_cumulative(), hurdle_cumulative())),
  sr_acat_th = list(f = bf(y ~ x, theta1 ~ z),
                    fam = mixture(sratio(), acat())),
  gr_cum2 = list(f = bf(y | thres(gr = gr) ~ x),
                 fam = mixture(cumulative(), cumulative())),
  # the third warning's guard has no family condition: components of
  # different families at shared thresholds are different distributions
  mu_cum_sr_theta = list(f = bf(ysh ~ 1, theta1 ~ z),
                         fam = mixture(cumulative(), sratio(),
                                       order = "mu")),
  # the hurdle warning's own remedy: priors on hu1 and hu2
  hurdle_prior = list(f = bf(yh ~ x),
                      fam = mixture(hurdle_cumulative(), hurdle_cumulative()),
                      prior = c(set_prior("beta(4, 16)", class = "hu1"),
                                set_prior("beta(8, 12)", class = "hu2"))),
  # the same shared-threshold case with two links of one family, and the
  # true-positive control: one family, one link, theta1 ~ z
  mu_cum_links_theta = list(f = bf(ysh ~ 1, theta1 ~ z),
                            fam = mixture(cumulative(), cumulative("probit"),
                                          order = "mu")),
  mu_same_theta = list(f = bf(ysh ~ 1, theta1 ~ z),
                       fam = mixture(cumulative(), cumulative(),
                                     order = "mu")),
  plain_cum = list(f = bf(y ~ x), fam = cumulative()),
  plain_probit = list(f = bf(ypr ~ x), fam = cumulative("probit"))
)
for (s in seq_len(R)) {
  d <- gen(s)
  warns <- character(0)
  err <- NULL
  fit <- tryCatch(withCallingHandlers(
    frm(spec$f, family = spec$fam, data = d, prior = spec$prior),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      err <<- conditionMessage(e)
      NULL
    })
  if (is.null(fit)) {
    cat(sprintf("REP model=%s seed=%d ERROR %s\n", model, s,
                substr(err, 1, 160)))
    next
  }
  gm <- max(abs(fit$obj$gr(fit$opt$par)))
  se <- tryCatch(withCallingHandlers(
    fixef(fit)[, "Est.Error"],
    warning = function(w) {
      warns <<- c(warns, paste("[se]", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }), error = function(e) NA_real_)
  flat <- grepl("flat|not identified|no usable standard error|one distribution|not in the likelihood",
                warns)
  cat(sprintf("REP model=%s seed=%d ll=%.6f conv=%d gradfinite=%s maxgrad=%.3g se_all_finite=%s maxse=%.3g nwarn=%d nflat=%d\n",
              model, s, as.numeric(logLik(fit)), fit$opt$convergence,
              is.finite(gm), gm, all(is.finite(se)), max(se), length(warns),
              sum(flat)))
  for (w in unique(warns)) cat("   WARN", substr(w, 1, 200), "\n")
}
