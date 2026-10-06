# Reviewer, re-check: defect 8, `ls ~ mo(income) * age` with NaN
# standard errors. Own construction: seeds 1 to 40 of brms_monotonic's
# data code; per seed the optimizer code, every warning AND message,
# whether any standard error is finite, the gradient at the optimum,
# and whether a second optimizer reaches the same log-likelihood (so a
# NaN standard error is not an unconverged fit in disguise).
#
#   Rscript dev/vigport-rev2-mo.R [nseeds]
args <- commandArgs(trailingOnly = TRUE)
ns <- if (length(args)) as.integer(args[1]) else 40L
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
quiet_fit <- function(...) {
  cond <- character()
  f <- withCallingHandlers(frm(...), warning = function(w) {
    cond <<- c(cond, paste("W:", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cond <<- c(cond, paste("M:", conditionMessage(m)))
    invokeRestart("muffleMessage")
  })
  list(f = f, cond = cond)
}
rows <- list()
for (s in seq_len(ns)) {
  d <- mk(s)
  r <- quiet_fit(ls ~ mo(income) * age, data = d)
  f <- r$f
  se <- fixef(f)[, "Est.Error"]
  r2 <- quiet_fit(ls ~ mo(income) * age, data = d,
                  control = frmtmb_control(optimizer = "optim"))
  g <- tryCatch(max(abs(f$obj$gr(f$opt$par))), error = function(e) NA)
  # a centered age: the same model, reparameterized
  d$agec <- d$age - mean(d$age)
  r3 <- quiet_fit(ls ~ mo(income) * agec, data = d)
  rows[[s]] <- data.frame(
    seed = s, code = f$opt$convergence, n_cond = length(r$cond),
    finite_se = sum(is.finite(se)), n_se = length(se),
    maxgrad = g, ll = as.numeric(logLik(f)),
    ll_optim = as.numeric(logLik(r2$f)),
    finite_se_centered = sum(is.finite(fixef(r3$f)[, "Est.Error"])),
    ll_centered = as.numeric(logLik(r3$f)),
    cond = substr(paste(r$cond, collapse = " | "), 1, 70))
}
X <- do.call(rbind, rows)
X$all_nan <- X$finite_se == 0
print(X[, c("seed", "code", "n_cond", "finite_se", "n_se", "maxgrad",
            "ll", "ll_optim", "finite_se_centered", "ll_centered")],
      row.names = FALSE, digits = 6)
cat(sprintf("\nseeds %d; all SEs non-finite on %d; of those code 0 and no",
            ns, sum(X$all_nan)))
cat(sprintf(" warning or message: %d\n",
            sum(X$all_nan & X$code == 0 & X$n_cond == 0)))
cat("some but not all SEs non-finite:",
    sum(X$finite_se > 0 & X$finite_se < X$n_se), "\n")
cat("max |logLik(nlminb) - logLik(optim)| on all-NaN seeds:",
    format(max(abs(X$ll - X$ll_optim)[X$all_nan])), "\n")
cat("centered age: all-NaN seeds with finite SEs after centering:",
    sum(X$all_nan & X$finite_se_centered == X$n_se), "of", sum(X$all_nan),
    "; max |logLik diff|", format(max(abs(X$ll - X$ll_centered))), "\n")
s0 <- X$seed[X$all_nan & X$code == 0 & X$n_cond == 0][1]
if (!is.na(s0)) {
  cat("\nsummary() at seed", s0, "(code 0, no condition):\n")
  print(summary(quiet_fit(ls ~ mo(income) * age, data = mk(s0))$f))
}
