# Lane splinecurve, round 2: what the ROUND-1 build answers on the two
# blockers, as numbers, before the round-2 fix is installed.
#   Rscript dev/splinecurve-r2-before.R > dev/splinecurve-r2-before-numbers.log
# Fixture: test-new-levels.R's own (set.seed(4), 12 subjects, 1200 rows).
lib <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline", format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "; formals of sp_spec:",
    paste(names(formals(frmtmb.spline:::sp_spec)), collapse = ", "), "\n\n")
set.seed(4)
n_sub <- 12; n_rep <- 4; n_t <- 25
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(subject = factor(sub),
                t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h <- rnorm(n_sub, 1, 0.12)
sv <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h[sub], sv[sub]) + rnorm(nrow(d), 0, 0.06)
tt <- seq(0.05, 0.95, length.out = 30)
gu <- data.frame(t = tt, subject = factor("new",
                                          levels = c(levels(d$subject),
                                                     "new")))
show <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  cat("--", label, ":", if (inherits(r, "error")) {
    paste("REFUSED", substr(conditionMessage(r), 1, 90))
  } else "ANSWERED", "\n")
  if (!inherits(r, "error")) r else NULL
}

## Blocker 1
for (form in list(v ~ s(t, k = 10) + (1 | subject),
                  v ~ s(t, k = 10) + (1 + t | subject))) {
  cat("====", deparse(form), "====\n")
  fit <- suppressWarnings(frm(bf(form), family = gaussian(), data = d))
  cv <- show("frm_curve(simultaneous = TRUE), unseen, NULL, TRUE",
             frm_curve(fit, newdata = gu, re_formula = NULL,
                       allow_new_levels = TRUE, nsim = 20000, seed = 1))
  if (!is.null(cv)) {
    cat("   crit_sim", format(cv$.crit_sim[1], digits = 5),
        "against qnorm(0.975)", format(qnorm(0.975), digits = 5), "\n")
  }
  show("frm_curve_deriv, unseen, NULL, TRUE",
       frm_curve_deriv(fit, var = "t", newdata = gu, re_formula = NULL,
                       allow_new_levels = TRUE, simultaneous = FALSE))
  pk <- show("frm_curve_feature(maximum), unseen, NULL, TRUE",
             frm_curve_feature(fit, var = "t", type = "maximum",
                               newdata = gu, re_formula = NULL,
                               allow_new_levels = TRUE))
  if (!is.null(pk)) {
    at <- frm_curve(fit, newdata = transform(gu[1, ], t = pk$.estimate),
                    re_formula = NULL, allow_new_levels = TRUE,
                    simultaneous = FALSE)
    cat("   .value_se", format(pk$.value_se, digits = 5),
        " frm_curve .se at the same t", format(at$.se, digits = 5), "\n")
  }
}

## Blocker 2
cat("==== stored-curve rule ====\n")
fit <- suppressWarnings(frm(bf(v ~ s(t, k = 10) +
                                 s(t, subject, bs = "fs", k = 5)),
                            family = gaussian(), data = d))
seen <- data.frame(t = tt, subject = factor("1", levels = levels(d$subject)))
st <- frm_curve(fit, newdata = seen, re_formula = NA,
                allow_new_levels = TRUE, simultaneous = FALSE)
show("frm_curve_deriv(stored TRUE, unseen grid, EXPLICIT FALSE)",
     frm_curve_deriv(st, var = "t", newdata = gu, allow_new_levels = FALSE,
                     simultaneous = FALSE))
show("frm_curve_feature(stored TRUE, unseen grid, EXPLICIT FALSE)",
     frm_curve_feature(st, var = "t", newdata = gu,
                       allow_new_levels = FALSE))
