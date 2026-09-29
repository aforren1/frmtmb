# Reviewer, lane splinecurve, re-check: the stored-curve rule, as the
# full matrix of stored value x caller's value, for deriv and feature.
#   Rscript splinecurve-rev-09-stored.R > splinecurve-rev-log/r2-09-stored.txt
# Fixture: test-new-levels.R's own (set.seed(4), 12 subjects, 1200 rows).
LIB <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
set.seed(4)
n_sub <- 12; n_rep <- 4; n_t <- 25
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(subject = factor(sub),
                t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h <- rnorm(n_sub, 1, 0.12)
sv <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h[sub], sv[sub]) + rnorm(nrow(d), 0, 0.06)
fit <- suppressWarnings(frm(bf(v ~ s(t, k = 10) + s(t, subject, bs = "fs",
                                                   k = 5)),
                            family = gaussian(), data = d))
tt <- seq(0.05, 0.95, length.out = 30)
lev <- c(levels(d$subject), "population")
unseen <- data.frame(t = tt, subject = factor("population", levels = lev))
seen <- data.frame(t = tt, subject = factor("1", levels = levels(d$subject)))
res <- function(expr) tryCatch(expr, error = function(e) e)
say <- function(label, r) {
  cat(sprintf("%-58s %s\n", label, if (inherits(r, "error")) {
    paste0("REFUSED [", class(r)[2], "]")
  } else paste("ANSWERED, rows", nrow(r))))
}
stored <- list(
  `TRUE` = frm_curve(fit, newdata = seen, allow_new_levels = TRUE,
                     simultaneous = FALSE),
  `FALSE` = frm_curve(fit, newdata = seen, simultaneous = FALSE))
ref_d <- frm_curve_deriv(fit, var = "t", newdata = unseen,
                         allow_new_levels = TRUE, simultaneous = FALSE)
for (st in names(stored)) {
  cv <- stored[[st]]
  for (cl in c("missing", "FALSE", "TRUE")) {
    dd <- if (cl == "missing") {
      res(frm_curve_deriv(cv, var = "t", newdata = unseen,
                          simultaneous = FALSE))
    } else {
      res(frm_curve_deriv(cv, var = "t", newdata = unseen,
                          allow_new_levels = as.logical(cl),
                          simultaneous = FALSE))
    }
    ff <- if (cl == "missing") {
      res(frm_curve_feature(cv, var = "t", type = "maximum",
                            newdata = unseen))
    } else {
      res(frm_curve_feature(cv, var = "t", type = "maximum",
                            newdata = unseen,
                            allow_new_levels = as.logical(cl)))
    }
    say(paste0("deriv:   stored ", st, ", caller ", cl, ", unseen grid"), dd)
    if (!inherits(dd, "error")) {
      cat("   identical() to deriv on the fit with TRUE:",
          identical(dd$.estimate, ref_d$.estimate) &&
            identical(dd$.se, ref_d$.se), "\n")
    }
    say(paste0("feature: stored ", st, ", caller ", cl, ", unseen grid"), ff)
  }
  # the stored grid itself (seen), caller FALSE: must still answer
  say(paste0("deriv:   stored ", st, ", caller FALSE, stored seen grid"),
      res(frm_curve_deriv(cv, var = "t", allow_new_levels = FALSE,
                          simultaneous = FALSE)))
}
# a stored curve AT the unseen level, re-read with the caller's FALSE
cu <- frm_curve(fit, newdata = unseen, allow_new_levels = TRUE,
                simultaneous = FALSE)
say("deriv:   stored TRUE at unseen, caller missing, stored grid",
    res(frm_curve_deriv(cu, var = "t", simultaneous = FALSE)))
say("deriv:   stored TRUE at unseen, caller FALSE, stored grid",
    res(frm_curve_deriv(cu, var = "t", allow_new_levels = FALSE,
                        simultaneous = FALSE)))
