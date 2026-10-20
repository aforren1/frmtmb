# Review of lane adhess, item 1: no-RE fits of frmtmb.learn, frmtmb.ode
# and frmtmb.latent against a Richardson Jacobian of the AD gradient.
#   REV_ARM=lane|base Rscript dev/adhess-rev-ext.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
suppressMessages({
  library(frmtmb.learn)
  library(frmtmb.ode)
  library(frmtmb.latent)
})
ext <- "C:/Users/adf44/source/r/frmtmb-wt-adhess/extensions/"
defs <- function(f) {
  for (e in parse(paste0(ext, f))) {
    if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, globalenv())
  }
}
defs("frmtmb.learn/tests/testthat/test-counterfactual.R")
defs("frmtmb.latent/tests/testthat/test-hmm-starts.R")
cases <- list(
  learn = function() {
    d <- ct_design()
    frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1),
        family = bandit2arm_delta(subject = id, trial = trial), data = d)
  },
  ode = function() {
    set.seed(4321)
    ns <- 25L
    ts <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
    d <- data.frame(id = factor(rep(seq_len(ns), each = length(ts))),
                    time = rep(ts, ns), dose = 100)
    ka <- 1; ke <- 0.2
    d$conc <- 100 * ka / (10 * (ka - ke)) *
      (exp(-ke * d$time) - exp(-ka * d$time)) +
      stats::rnorm(nrow(d), 0, 0.2)
    frm(bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                          V = exp(lV)),
                             times = time, group = id, ncmt = 1,
                             depot = TRUE, init = list(depot = dose)),
           lka ~ 1, lke ~ 1, lV ~ 1, nl = TRUE) + gaussian(),
        data = d, start = list(beta = c(0, log(0.25), log(8))))
  },
  latent = function() {
    dd <- small_data(seed = 4501L)
    frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
        data = dd)
  })
for (nm in names(cases)) {
  r <- cap(cases[[nm]]())
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(nm, "ERROR", fit, "\n")
    next
  }
  at <- the_point(fit)
  V <- suppressWarnings(vcov(fit, full = TRUE))
  Hr <- rich_hess(fit$obj, at)
  e <- cov_err(V, solve(Hr))
  h <- tryCatch(ns$fit_outer_hessian(fit), error = function(e) NULL)
  cat(sprintf(paste0("%-7s np %d re %d | exact %s | vcov vs Richardson se ",
                     "%.2e cor %.2e | lost %s | warn %d\n"),
              nm, length(at), length(fit$obj$env$random),
              format(isTRUE(h$exact)), e[["se"]], e[["cor"]],
              paste(names(sdr_of(fit)$se_lost %||% character(0)),
                    collapse = ","), length(r$w)))
}
