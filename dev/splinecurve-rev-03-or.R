# Reviewer, lane splinecurve, claim 2: the stored-curve OR rule.
#   Rscript splinecurve-rev-03-or.R <lib> > splinecurve-rev-log/03-or-<lib>.txt
# Fixture: test-new-levels.R's own (set.seed(4), 12 subjects, 1200 rows),
# plus the same data with (1 | subject) in place of the fs term.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline", format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "\n\n")
has_flag <- "allow_new_levels" %in% names(formals(frm_curve))
res <- function(expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    paste0("REFUSED [", class(r)[2], "]: ", substr(conditionMessage(r), 1, 90))
  } else r
}
show <- function(label, r) {
  if (is.character(r)) cat(label, ":", r, "\n") else
    cat(label, ": ANSWERED, rows", nrow(r), "\n")
  invisible(r)
}

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
fit_re <- frm(bf(v ~ s(t, k = 10) + (1 | subject)), family = gaussian(),
              data = d)
tt <- seq(0.05, 0.95, length.out = 30)
lev <- c(levels(d$subject), "population")
unseen <- data.frame(t = tt, subject = factor("population", levels = lev))
seen <- data.frame(t = tt, subject = factor("1", levels = levels(d$subject)))

cat("-- before-arm check: the default refuses an unseen level of both",
    "models\n")
show("frm_curve(fs fit, unseen)", res(frm_curve(fit, newdata = unseen,
                                               simultaneous = FALSE)))
show("frm_curve(re fit, unseen, re_formula = NULL)",
     res(frm_curve(fit_re, newdata = unseen, re_formula = NULL,
                   simultaneous = FALSE)))
show("frm_curve(re fit, no subject column, re_formula = NULL)",
     res(frm_curve(fit_re, newdata = unseen[, "t", drop = FALSE],
                   re_formula = NULL, simultaneous = FALSE)))
if (!has_flag) quit(save = "no")

cat("\n-- a curve stored with TRUE on a grid of SEEN levels\n")
cs <- frm_curve(fit, newdata = seen, allow_new_levels = TRUE,
                simultaneous = FALSE)
cs0 <- frm_curve(fit, newdata = seen, simultaneous = FALSE)
cat("stored TRUE on a seen grid changes no number:",
    identical(cs$.estimate, cs0$.estimate) && identical(cs$.se, cs0$.se),
    "\n")
show("deriv(stored TRUE, newdata = unseen), caller default",
     res(frm_curve_deriv(cs, var = "t", newdata = unseen,
                         simultaneous = FALSE)))
show("deriv(stored TRUE, newdata = unseen), caller EXPLICIT FALSE",
     res(frm_curve_deriv(cs, var = "t", newdata = unseen,
                         allow_new_levels = FALSE, simultaneous = FALSE)))
show("deriv(FIT, newdata = unseen), default",
     res(frm_curve_deriv(fit, var = "t", newdata = unseen,
                         simultaneous = FALSE)))
cat("for comparison, re_formula on a stored curve: the caller's value is",
    "ignored\n")
cna <- frm_curve(fit_re, newdata = seen, simultaneous = FALSE)
dn <- frm_curve_deriv(cna, var = "t", re_formula = NULL, simultaneous = FALSE)
cat("  stored NA, caller NULL -> used:",
    deparse(attr(dn, "spec")$re_formula), "\n")

cat("\n-- the same chain on (1 | subject), where an unseen level carries",
    "variance\n")
csr <- frm_curve(fit_re, newdata = seen, re_formula = NULL,
                 allow_new_levels = TRUE, simultaneous = FALSE)
fr <- res(frm_curve_feature(csr, var = "t", type = "crossing", at = 0.5,
                            newdata = unseen))
show("feature(stored TRUE at a seen level, newdata = unseen), default", fr)
if (!is.character(fr)) {
  at_root <- frm_curve(fit_re, newdata = transform(unseen[seq_len(nrow(fr)), ],
                                                   t = fr$.estimate),
                       re_formula = NULL, allow_new_levels = TRUE,
                       simultaneous = FALSE)
  cat("  .value_se", format(fr$.value_se, digits = 5),
      " frm_curve .se at the root", format(at_root$.se, digits = 5), "\n")
}
