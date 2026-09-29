# Reviewer, lane splinecurve, re-check after punch round 1: every route
# that gave a wrong number at an unseen bar-term level refuses, the
# pointwise band still answers with the right se, and nothing refuses
# where extra_var is 0.
#   Rscript splinecurve-rev-08-recheck.R > splinecurve-rev-log/r2-08-recheck.txt
# Data: the vignette simulation, set.seed(4), as in rev-02.
LIB <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline", format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "\n\n")
set.seed(4)
n_sub <- 20; n_rep <- 12; n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(subject = factor(sub),
                t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
lev_u <- c(levels(d$subject), "new")
grid <- data.frame(t = seq(0.05, 0.95, length.out = 60),
                   subject = factor("new", levels = lev_u))
# a window on the far tail, where the curve does not cross 0.5
tail_grid <- data.frame(t = seq(0.9, 1, length.out = 20),
                        subject = factor("new", levels = lev_u))

res <- function(expr) tryCatch(expr, error = function(e) e)
say <- function(label, r) {
  if (inherits(r, "error")) {
    cat(sprintf("  %-46s REFUSED: %s\n", label,
                substr(conditionMessage(r), 1, 70)))
  } else {
    cat(sprintf("  %-46s ANSWERED, %d rows\n", label, nrow(r)))
  }
  invisible(r)
}
routes <- function(fit, nd, rf, anl, tag) {
  cat("--", tag, "re_formula =", deparse(rf), "allow_new_levels =", anl,
      "\n")
  lb <- res(frm_lp_basis(fit, newdata = nd, re_formula = rf,
                         allow_new_levels = anl))
  if (!inherits(lb, "error")) {
    cat("  extra_var range", format(range(lb$extra_var), digits = 5), "\n")
  }
  out <- list(
    pw = say("frm_curve(simultaneous = FALSE)",
             res(frm_curve(fit, newdata = nd, re_formula = rf,
                           allow_new_levels = anl, simultaneous = FALSE))),
    sim = say("frm_curve(simultaneous = TRUE)",
              res(frm_curve(fit, newdata = nd, re_formula = rf,
                            allow_new_levels = anl, nsim = 5000, seed = 1))),
    sim_tr = say("frm_curve(simultaneous = TRUE, transform = TRUE)",
                 res(frm_curve(fit, newdata = nd, re_formula = rf,
                               allow_new_levels = anl, nsim = 5000, seed = 1,
                               transform = TRUE))),
    d1 = say("frm_curve_deriv(order = 1)",
             res(frm_curve_deriv(fit, var = "t", newdata = nd,
                                 re_formula = rf, allow_new_levels = anl,
                                 simultaneous = FALSE))),
    d2 = say("frm_curve_deriv(order = 2)",
             res(frm_curve_deriv(fit, var = "t", order = 2, newdata = nd,
                                 re_formula = rf, allow_new_levels = anl,
                                 simultaneous = FALSE))),
    fx = say("frm_curve_feature(maximum)",
             res(frm_curve_feature(fit, var = "t", type = "maximum",
                                   newdata = nd, re_formula = rf,
                                   allow_new_levels = anl))),
    fc = say("frm_curve_feature(crossing, at = 0.5)",
             res(frm_curve_feature(fit, var = "t", type = "crossing",
                                   at = 0.5, newdata = nd, re_formula = rf,
                                   allow_new_levels = anl))),
    ft = say("feature(crossing 0.5) on a ROOTLESS tail window",
             res(frm_curve_feature(fit, var = "t", type = "crossing",
                                   at = 0.5, newdata = tail_grid,
                                   re_formula = rf,
                                   allow_new_levels = anl))))
  if (!inherits(lb, "error") && !inherits(out$pw, "error")) {
    C <- as.matrix(lb$A)
    se_right <- sqrt(diag(C %*% lb$V %*% t(C)) + lb$extra_var)
    ref <- frm_linpred(fit, newdata = nd, re_formula = rf,
                       allow_new_levels = anl, se.fit = TRUE)
    cat("  pointwise .se vs sqrt(diag(A V A') + extra_var): max rel",
        format(max(abs(out$pw$.se / se_right - 1)), digits = 3),
        "; vs frm_linpred(se.fit): max rel",
        format(max(abs(out$pw$.se / as.numeric(ref$se.fit) - 1)),
               digits = 3), "\n")
  }
  invisible(out)
}

for (form in list(v ~ s(t, k = 12) + (1 | subject),
                  v ~ s(t, k = 12) + (1 + t | subject))) {
  cat("====", deparse(form), "====\n")
  fit <- frm(bf(form), family = gaussian(), data = d)
  routes(fit, grid, NULL, TRUE, "unseen level")
  routes(fit, grid[, "t", drop = FALSE], NULL, TRUE, "grouping column absent")
  a <- routes(fit, grid, NA, TRUE, "unseen level")
  # at NA the grouping term is dropped, so TRUE must change no number
  # against a grid that names a SEEN subject
  seen <- transform(grid, subject = factor("3", levels = levels(d$subject)))
  b <- list(
    pw = frm_curve(fit, newdata = seen, re_formula = NA,
                   simultaneous = FALSE),
    sim = frm_curve(fit, newdata = seen, re_formula = NA, nsim = 5000,
                    seed = 1),
    d1 = frm_curve_deriv(fit, var = "t", newdata = seen, re_formula = NA,
                         simultaneous = FALSE),
    fc = frm_curve_feature(fit, var = "t", type = "crossing", at = 0.5,
                           newdata = seen, re_formula = NA))
  for (k in names(b)) {
    if (inherits(a[[k]], "error")) next
    cat("  NA, unseen+TRUE vs seen+FALSE,", k, ": .estimate and .se",
        "identical():", identical(a[[k]]$.estimate, b[[k]]$.estimate) &&
          identical(a[[k]]$.se, b[[k]]$.se), "\n")
  }
  cat("\n")
}

cat("==== vignette fs model, population curve at an unseen subject ====\n")
fit_fs <- frm(bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
              family = gaussian(), data = d)
pop <- factor("population", levels = c(levels(d$subject), "population"))
g80 <- data.frame(t = seq(0, 1, length.out = 80), subject = pop)
g2 <- data.frame(t = seq(0.05, 0.95, length.out = 40), subject = pop)
cv <- say("frm_curve(nsim = 20000, seed = 1)",
          res(frm_curve(fit_fs, newdata = g80, re_formula = NA,
                        allow_new_levels = TRUE, nsim = 20000, seed = 1)))
if (!inherits(cv, "error")) {
  cat("  crit", format(cv$.crit_sim[1], digits = 10), " mcse",
      format(attr(cv, "check")$crit_mcse, digits = 7), "\n")
}
d1 <- say("frm_curve_deriv(nsim = 20000, seed = 2)",
          res(frm_curve_deriv(fit_fs, var = "t", newdata = g2,
                              re_formula = NA, allow_new_levels = TRUE,
                              nsim = 20000, seed = 2)))
pk <- say("frm_curve_feature(maximum)",
          res(frm_curve_feature(fit_fs, var = "t", type = "maximum",
                                newdata = g2, re_formula = NA,
                                allow_new_levels = TRUE)))
if (!inherits(pk, "error")) {
  cat("  peak", format(pk$.estimate, digits = 7), " se",
      format(pk$.se, digits = 7), "\n")
}
say("frm_curve_feature(crossing, at = 0.2)",
    res(frm_curve_feature(fit_fs, var = "t", type = "crossing", at = 0.2,
                          newdata = g2, re_formula = NA,
                          allow_new_levels = TRUE)))
say("frm_curve(simultaneous = FALSE)",
    res(frm_curve(fit_fs, newdata = g80, re_formula = NA,
                  allow_new_levels = TRUE, simultaneous = FALSE)))

cat("\n==== exact gp(), TRUE, no unseen level (the other extra_var) ====\n")
set.seed(8)
dgp <- data.frame(x = sort(runif(60, 0, 5)))
dgp$y <- sin(dgp$x) + rnorm(60, 0, 0.3)
fgp <- frm(bf(y ~ gp(x)), family = gaussian(), data = dgp)
gx <- data.frame(x = dgp$x[-1] - diff(dgp$x) / 2)
cat("  extra_var range",
    format(range(frm_lp_basis(fgp, newdata = gx)$extra_var), digits = 3),
    "\n")
for (anl in c(FALSE, TRUE)) {
  cat("  allow_new_levels =", anl, "\n")
  say("  frm_curve(simultaneous = TRUE)",
      res(frm_curve(fgp, newdata = gx, allow_new_levels = anl, nsim = 2000,
                    seed = 1)))
  say("  frm_curve_deriv()",
      res(frm_curve_deriv(fgp, var = "x", newdata = gx,
                          allow_new_levels = anl, simultaneous = FALSE)))
}
