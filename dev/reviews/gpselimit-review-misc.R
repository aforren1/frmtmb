# Review of lane gpselimit, checks 5 and 7.
#  K: check_se = "warning" (default), "ignore", "stop" on a singular
#     (1 | g) fit and on a white gp() fit: what frm() and summary() say.
#  X: the lane's MEDIUM on exact gp() (dev/gpselimit-exact2.R's design,
#     rep 3 seeds 1 to 4) in a fresh process: summary()'s lscale row, the
#     lost reasons and the conditions.
# Usage: Rscript gpselimit-review-misc.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-review-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]
rv_libs(arm)
cat("check_se default:", format(frmtmb_control()$check_se), "\n")
set.seed(31)
dg <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
dg$y <- 1 + 0.5 * dg$x + rnorm(60)
s3 <- rv_sweep_row(3)
for (cs in c("warning", "ignore", "stop")) {
  for (case in c("singular", "white")) {
    q <- if (case == "singular") {
      rv_fit(frm(y ~ x + (1 | g), data = dg,
                 control = frmtmb_control(check_se = cs)))
    } else {
      rv_fit(frm(bf(y ~ gp(x, k = 10)), data = s3$d, family = gaussian(),
                 control = frmtmb_control(check_se = cs)))
    }
    if (inherits(q$fit, "error")) {
      cat(sprintf("K %s check_se=%s %s | frm ERROR: %s\n", arm, cs, case,
                  substr(conditionMessage(q$fit), 1, 90)))
      next
    }
    sm <- rv_fit(capture.output(print(summary(q$fit))))
    cat(sprintf("K %s check_se=%s %s | frm conds [%s] | summary conds [%s]%s\n",
                arm, cs, case, rv_cstr(q$conds), rv_cstr(sm$conds),
                if (inherits(sm$fit, "error")) {
                  paste0(" summary ERROR: ",
                         substr(conditionMessage(sm$fit), 1, 80))
                } else ""))
  }
}
for (seed in 1:4) {
  rep <- 3
  set.seed(7700 + 10 * rep + seed)
  np <- 40
  xs <- sort(runif(np, 0, 10))
  u <- rnorm(np)
  x <- rep(xs, each = rep)
  y <- rep(u, each = rep) + 0.3 * rnorm(np * rep)
  q <- rv_fit(frm(bf(y ~ gp(x)), data = data.frame(x = x, y = y),
                  family = gaussian()))
  fit <- q$fit
  sm <- rv_fit(capture.output(print(summary(fit))))
  row <- grep("lscale", sm$fit, value = TRUE)
  lost <- fit$cache$sdr$se_lost
  cat(sprintf("X %s rep3 seed %d | ll %.6f | frm conds [%s] | summary conds [%s] | lost [%s] | %s\n",
              arm, seed, as.numeric(logLik(fit)), rv_cstr(q$conds),
              rv_cstr(sm$conds),
              paste(paste0(names(lost), "=", lost), collapse = ","),
              paste(trimws(row), collapse = " / ")))
}
