# Reviewer: tier 2 on a model WITH random effects, where the Hessian is
# optimHess()'s finite differences (no exact obj$he()). A saturated mo()
# simplex coordinate has a row of order 1e-20 whose finite-difference
# noise, scaled to unit diagonal, could reach order 1. Reference: the
# same optimHess() matrix inverted with the saturated coordinates (weight
# below 1e-8) removed.
#   Rscript dev/nanse-rev-mo-re.R base|lane [seeds]
args <- commandArgs(TRUE)
arm <- args[1]
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:60
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d$grp <- factor(rep(1:10, 10))
  d$ls <- d$ls + rnorm(10, 0, 2)[d$grp]
  d
}
tot <- c(n = 0, allfinite = 0, tier2 = 0, worst = 0)
for (s in seeds) {
  d <- mk(s)
  w <- character()
  fit <- withCallingHandlers(frm(ls ~ mo(income) * age + (1 | grp), data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (fit$opt$convergence != 0) next
  nm <- ns$outer_par_names(fit)
  p <- fit$opt$par
  z <- grepl("^zeta", nm)
  # simplex weights from the softmax coordinates, per simplex
  sat <- rep(FALSE, length(p))
  for (pre in unique(sub("_[0-9]+$", "", nm[z]))) {
    i <- which(sub("_[0-9]+$", "", nm) == pre & z)
    e <- exp(c(0, p[i]) - max(c(0, p[i])))
    wgt <- e / sum(e)
    sat[i] <- wgt[-1] < 1e-8
  }
  if (!any(sat)) next
  sdr <- ns$sdr_of(fit)
  raw <- RTMB::sdreport(fit$obj)$cov.fixed
  plain_ok <- all(is.finite(sqrt(diag(raw))))
  se <- suppressWarnings(sqrt(diag(sdr$cov.fixed)))
  H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
  keep <- !sat
  ref <- sqrt(diag(solve(H[keep, keep])))
  b <- nm %in% c("(Intercept)", "age", "moincome", "moincome:age",
                 "sigma_(Intercept)", "theta_1")
  r <- max(abs(se[b & keep] / ref[b[keep]] - 1))
  tot["n"] <- tot["n"] + 1
  if (!plain_ok && all(is.finite(se))) tot["tier2"] <- tot["tier2"] + 1
  if (all(is.finite(se[b]))) tot["allfinite"] <- tot["allfinite"] + 1
  if (is.finite(r)) tot["worst"] <- max(tot["worst"], r)
  cat(sprintf(paste0("seed %3d: saturated %s; plain SEs finite %s; ",
                     "reported coef SEs finite %s; max |SE/ref - 1| %.3g; ",
                     "warnings %d%s\n"),
              s, paste(nm[sat], collapse = ","), plain_ok,
              all(is.finite(se[b])), r, length(w),
              if (length(sdr$se_lost)) paste0(" lost ",
                paste(names(sdr$se_lost), collapse = ",")) else ""))
}
print(tot)
