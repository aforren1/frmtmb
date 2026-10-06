# Punch round 1, B2: the fits of the test suites where the degenerate
# warning fired (dev/ordmix-p1-check/suite-log/), refitted on the lane
# build, each with the review's label (any |estimate| > 30 or a
# non-finite standard error) and the check's own numbers. Was any of
# them a sound fit?
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < stats::plogis(-1.2 + 0.5 * d$z), 0L, d$y)
  d
}
alt_data <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < 0.2, 0L, d$y)
  d
}
gen19 <- function() {
  set.seed(19)
  x <- rnorm(300)
  cls <- rbinom(300, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(300)
  data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
}
st <- "student_t(3, 0, 2.5)"
cases <- list(
  simulate = function() frm(bf(y | thres(gr = g) ~ x),
                            family = mixture(cumulative(), sratio()),
                            data = omx_data(20261011, n = 300)),
  refusals_cs = function() frm(bf(y ~ x, mu2 ~ cs(x)),
                               family = mixture(cumulative(), sratio()),
                               data = omx_data(20261012)),
  b2_seed19 = function() frm(bf(y ~ x), family = mixture(cumulative(),
                                                         cumulative()),
                             data = gen19()),
  thres4 = function() frm(bf(y | thres(4) ~ x),
                          family = mixture(cumulative(), sratio()),
                          data = omx_data(20261016)),
  brms_cs = function() frm(bf(y ~ cs(x)), family = mixture(sratio(), acat()),
                           data = alt_data(20261030, 300)),
  draws_none = function() frm(bf(y ~ x, disc1 ~ 0 + z),
    family = mixture(cumulative(), sratio()), data = alt_data(71, 200),
    prior = set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
      set_prior("normal(0, 2)", class = "b", dpar = "mu2") +
      set_prior(st, class = "Intercept", dpar = "mu1") +
      set_prior(st, class = "Intercept", dpar = "mu2")))
for (nm in names(cases)) {
  w <- character(0)
  fit <- withCallingHandlers(cases[[nm]](), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  fx <- suppressWarnings(fixef(fit))
  dg <- frmtmb:::mixture_ord_degeneracy(fit, names(fit$spec$responses)[1])
  cat(sprintf(paste0("FIRE case=%s label=%s max|est|=%.4g se_finite=%s ",
                     "reach=%s collapsed=%s degenerate_warning=%s\n"),
              nm, any(abs(fx[, "Estimate"]) > 30) ||
                any(!is.finite(fx[, "Est.Error"])),
              max(abs(fx[, "Estimate"])), all(is.finite(fx[, "Est.Error"])),
              paste(signif(dg$reach, 4), collapse = ","),
              paste(dg$collapsed, collapse = ","),
              any(grepl("degenerate boundary", w, fixed = TRUE))))
}
