# Reviewer, punch round 1 (B2 gap): frmtmb.spline's frm_curve() along a
# lost direction. y ~ a + b with a ~ 1 + x and b ~ 1 + x: only a + b is
# identified, so a curve of `a` alone over x has no standard error,
# while the curve of mu does. Compares frm_curve()'s band with
# frm_lp_basis()'s se_nonest flag and with fitted(dpar = "a").
#   Rscript dev/nanse-rev2-curve.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spline)})
cat("arm", arm, "frmtmb", find.package("frmtmb"), "spline",
    find.package("frmtmb.spline"), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "err")
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
set.seed(955)
n <- 60
dd <- data.frame(x = rnorm(n))
dd$y <- 3 + 0.5 * dd$x + rnorm(n, 0, 0.4)
fit <- cap(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE), data = dd))
cat("frm() warnings:", length(fit$w), "\n")
for (x in fit$w) cat("   ", substr(x, 1, 150), "\n")
fit <- fit$v
nd <- data.frame(x = seq(-2, 2, length.out = 9))
for (dp in c("mu", "a")) {
  dpa <- if (dp == "mu") NULL else dp
  lb <- cap(frm_lp_basis(fit, newdata = nd, dpar = dpa))
  nonest <- if (inherits(lb$v, "err")) NA else sum(lb$v$se_nonest)
  for (sim in c(FALSE, TRUE)) {
    cv <- cap(frm_curve(fit, newdata = nd, dpar = dpa,
                        simultaneous = sim, nsim = 2000, seed = 1))
    if (inherits(cv$v, "err")) {
      cat(sprintf("dpar %-2s simultaneous %-5s: ERROR %s\n", dp, sim,
                  substr(cv$v, 1, 600)))
      next
    }
    se <- cv$v$.se
    cat(sprintf(paste0("dpar %-2s simultaneous %-5s: se_nonest rows %s of ",
                       "%d; frm_curve .se %s; finite %d of %d; ",
                       "warnings %d%s\n"),
                dp, sim, nonest, nrow(nd),
                paste(signif(range(se), 3), collapse = " to "),
                sum(is.finite(se)), length(se), length(cv$w),
                if (length(cv$w)) paste0(" [", substr(cv$w[1], 1, 60), "]")
                else ""))
  }
  # the derivative of `a` along x is a_x, which the data do not determine
  for (sim in c(FALSE, TRUE)) {
    dv <- cap(frm_curve_deriv(fit, var = "x", newdata = nd, dpar = dpa,
                              simultaneous = sim, nsim = 2000, seed = 1))
    if (inherits(dv$v, "err")) {
      cat(sprintf("dpar %-2s deriv simultaneous %-5s: ERROR %s\n", dp, sim,
                  substr(dv$v, 1, 300)))
    } else {
      cat(sprintf(paste0("dpar %-2s deriv simultaneous %-5s: .se %s; ",
                         "finite %d of %d; warnings %d%s\n"), dp, sim,
                  paste(signif(range(dv$v$.se), 3), collapse = " to "),
                  sum(is.finite(dv$v$.se)), nrow(dv$v), length(dv$w),
                  if (length(dv$w)) paste0(" [", substr(dv$w[1], 1, 60), "]")
                  else ""))
    }
  }
  ft <- cap(fitted(fit, newdata = nd, dpar = dpa))
  cat(sprintf("dpar %-2s fitted() Est.Error finite %d of %d, warnings %d\n",
              dp, sum(is.finite(ft$v[, 2])), nrow(nd), length(ft$w)))
}
