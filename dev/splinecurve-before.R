# Lane splinecurve: the BEHAVIORAL failure on the unfixed package.
#
# Run against the read-only release library, which holds frmtmb 0.65.0
# and frmtmb.spline 0.8.1 as it stood before this lane:
#   Rscript dev/splinecurve-before.R > dev/splinecurve-before.log
# Data: the curve-inference vignette's own simulation, set.seed(4).
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
cat("frmtmb.spline", format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "\n")
cat("frm_curve() formals:", paste(names(formals(frm_curve)), collapse = ", "),
    "\n\n")

set.seed(4)
n_sub <- 20
n_rep <- 12
n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(
  subject = factor(sub),
  trial = rep(seq_len(n_sub * n_rep), each = n_t),
  t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
fit <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
           family = gaussian(), data = d)

try_it <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  cat("--", label, "\n")
  if (inherits(r, "error")) {
    cat("   ERROR [", paste(class(r), collapse = "/"), "]: ",
        conditionMessage(r), "\n\n", sep = "")
  } else {
    cat("   returned\n\n")
  }
}

# 1. the vignette's own call, as it stood
grid <- data.frame(t = seq(0, 1, length.out = 80))
try_it("vignette chunk `curve`: grid without subject, re_formula = NA",
       frm_curve(fit, newdata = grid, re_formula = NA, nsim = 20000,
                 seed = 1))
# 2. the documented population idiom, which 0.8.1 has no way to request
gu <- data.frame(t = grid$t,
                 subject = factor("population",
                                  levels = c(levels(d$subject),
                                             "population")))
try_it("grid at an UNSEEN subject, re_formula = NA (no flag to pass)",
       frm_curve(fit, newdata = gu, re_formula = NA, nsim = 20000, seed = 1))
try_it("grid at an UNSEEN subject, re_formula = NULL",
       frm_curve(fit, newdata = gu, re_formula = NULL, nsim = 20000,
                 seed = 1))
try_it("frm_curve_deriv at an UNSEEN subject",
       frm_curve_deriv(fit, var = "t", newdata = gu, re_formula = NA,
                       simultaneous = FALSE))
try_it("frm_curve_feature at an UNSEEN subject",
       frm_curve_feature(fit, var = "t", newdata = gu, re_formula = NA))
# 3. core itself answers the same question, so the gap is the extension's
p <- frm_linpred(fit, newdata = gu, re_formula = NA, allow_new_levels = TRUE)
cat("-- core frm_linpred(allow_new_levels = TRUE) at the unseen subject:",
    length(p), "values, first", format(p[1]), "\n")
