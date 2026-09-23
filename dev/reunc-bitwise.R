# Lane wt-reunc: what must NOT move. Run once per library and compare
# the two RDS files with dev/reunc-bitwise-cmp.R.
#
#   Rscript dev/reunc-bitwise.R <lib> <out.rds>
#
# Every fit is built in each mode the task names: ML, REML = TRUE and
# control(profile = TRUE), on a gaussian, a poisson and a distributional
# (sigma ~ x) model, each with a group-level term, plus two models the
# group-effect draw must leave alone because it has nothing to draw: no
# group-level term, and a population smooth only.
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
out <- a[2]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")

set.seed(20260922)
G <- 12; m <- 6; n <- G * m
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(n),
                t = runif(n))
u <- rnorm(G, 0, 0.6)
d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(n, 0, exp(-0.2 + 0.3 * d$x))
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x + u[d$g]))
d$ys <- sin(2 * pi * d$t) + rnorm(n, 0, 0.5)
nd <- data.frame(x = c(-1, 0, 1.5), t = c(0.2, 0.5, 0.8),
                 g = factor(c(1, 5, 9), levels = levels(d$g)))
ndnew <- data.frame(x = c(-1, 0, 1.5), g = factor(c("n1", "n2", "n1")))

models <- list(
  gauss = list(f = bf(y ~ x + (1 | g)), fam = gaussian()),
  pois = list(f = bf(cnt ~ x + (1 | g)), fam = poisson()),
  dist = list(f = bf(y ~ x + (1 | g), sigma ~ x), fam = gaussian()),
  nore = list(f = bf(y ~ x), fam = gaussian()),
  smooth = list(f = bf(ys ~ s(t, k = 6)), fam = gaussian())
)
modes <- list(
  ML = list(),
  REML = list(REML = TRUE),
  profile = list(control = frmtmb_control(profile = TRUE))
)
res <- list()
grab <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) {
    structure(conditionMessage(e), class = "grab_error")
  })
}
for (mn in names(models)) {
  for (md in names(modes)) {
    # poisson has no REML (REML integrates the gaussian mu coefficients)
    key <- paste(mn, md, sep = ".")
    mdl <- models[[mn]]
    fit <- grab(do.call(frm, c(list(mdl$f + mdl$fam, data = d),
                               modes[[md]])))
    if (inherits(fit, "grab_error")) {
      res[[key]] <- list(fit_error = unclass(fit))
      next
    }
    hasre <- mn %in% c("gauss", "pois", "dist")
    r <- list()
    r$vcov <- grab(vcov(fit))
    r$vcov_full <- grab(vcov(fit, full = TRUE))
    r$fixef <- grab(fixef(fit))
    r$summary_fixed <- grab(summary(fit)$fixed)
    r$summary_spec <- grab(summary(fit)$spec_pars)
    r$summary_random <- grab(summary(fit)$random)
    r$logLik <- grab(as.numeric(logLik(fit)))
    # the population level, which the task says must not move
    r$lp_na <- grab(frm_linpred(fit, re_formula = NA, se.fit = TRUE))
    r$lp_na_nd <- grab(frm_linpred(fit, newdata = nd, re_formula = NA,
                                   se.fit = TRUE))
    r$fitted_na <- grab(fitted(fit, re_formula = NA))
    r$fitted_na_nd <- grab(fitted(fit, newdata = nd, re_formula = NA))
    set.seed(11)
    r$predict_na <- grab(predict(fit, re_formula = NA, ndraws = 200))
    set.seed(12)
    r$predict_na_nd <- grab(predict(fit, newdata = nd, re_formula = NA,
                                    ndraws = 200, summary = FALSE))
    set.seed(13)
    r$predict_na_plug <- grab(predict(fit, newdata = nd, re_formula = NA,
                                      ndraws = 200,
                                      propagate_error = FALSE))
    # after predict(): the caller's stream must be left where it was
    r$stream_after <- stats::runif(1)
    # the conditional point values and the scalar standard errors: the
    # scalar delta method already read the joint covariance, so neither
    # the estimate nor its standard error has a reason to move
    r$lp_null <- grab(frm_linpred(fit, se.fit = TRUE))
    r$fitted_null <- grab(fitted(fit))
    r$fitted_null_nd <- grab(fitted(fit, newdata = nd))
    if (hasre) {
      set.seed(14)
      r$predict_new <- grab(predict(fit, newdata = ndnew, ndraws = 200,
                                    allow_new_levels = TRUE,
                                    summary = FALSE))
      r$fitted_new <- grab(fitted(fit, newdata = ndnew,
                                  allow_new_levels = TRUE))
      set.seed(15)
      # the conditional predict(), which the task DOES change; recorded
      # so the comparison can say it moved
      r$predict_null <- grab(predict(fit, newdata = nd, ndraws = 200))
      # and the caller's STREAM after each call that draws group
      # effects. A reviewer's rebuild found 9 differing stream positions
      # this battery never asked about, on rows whose predictions were
      # identical: the drawer took its seeds before predict() captured
      # the state it restores (punch round 1, M2)
      set.seed(21)
      invisible(grab(predict(fit, newdata = nd, ndraws = 200)))
      r$stream_after_cond <- stats::runif(1)
      set.seed(22)
      invisible(grab(predict(fit, newdata = ndnew, ndraws = 200,
                             allow_new_levels = TRUE)))
      r$stream_after_newlev <- stats::runif(1)
      set.seed(23)
      invisible(grab(predict(fit, ndraws = 100, summary = FALSE)))
      r$stream_after_insample <- stats::runif(1)
      set.seed(24)
      invisible(grab(fitted(fit, newdata = nd)))
      r$stream_after_fitted <- stats::runif(1)
    } else {
      set.seed(16)
      r$predict_null <- grab(predict(fit, newdata = nd, ndraws = 200))
      set.seed(17)
      r$predict_null_insample <- grab(predict(fit, ndraws = 100,
                                              summary = FALSE))
    }
    res[[key]] <- r
  }
}
saveRDS(res, out)
cat("wrote", out, "with", length(res), "fits\n")
