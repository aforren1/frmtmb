# Round 1b, item 2: the NaN-to-Inf trial mapping, checked on the lane
# build that has it.
#
#   Rscript dev/famlink-1b-nancheck.R > dev/famlink-1b-nancheck-log.txt
#
# Each design is fitted twice in ONE process on the same build:
#   old  nlminb driven on the raw objective through frmtmb_control's
#        custom-optimizer hook, which is the unmapped path the package
#        took before (same start, same control, scale 1)
#   new  the package's own nlminb path, which maps a NaN at a trial point
#        to +Inf
# and compared on NaN warnings, logLik and every outer estimate with
# identical(). Designs: the four of dev/famlink-punch-invgauss.R (seeds
# there), one whose optimum sits near eta = 0 (seed 20260917), and the
# broken custom family of part (c).
.libPaths(c("C:/Users/adf44/source/r/famlink-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

src <- parse("dev/famlink-punch-invgauss.R", keep.source = FALSE)
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(as.character(e[[2]]), "designs")) eval(e, globalenv())
}
designs$near_zero <- function() {
  # 1/mu^2 = eta: an intercept of 0.02 is a mean of about 7, and the
  # slope keeps some rows' eta within 0.01 of zero at the truth
  set.seed(20260917)
  n <- 400
  x <- runif(n, -1, 1)
  eta <- 0.02 + 0.0095 * x
  data.frame(y = RTMBdist::rinvgauss(n, 1 / sqrt(eta), 50), x = x)
}

raw_nlminb <- function(par, fn, gr, lower, upper, control) {
  r <- stats::nlminb(par, fn, gr, control = control, lower = lower,
                     upper = upper)
  r
}

fit_once <- function(d, family, arm) {
  warns <- character(0)
  ctl <- if (arm == "old") frmtmb_control(optimizer = raw_nlminb) else
    frmtmb_control()
  fit <- withCallingHandlers(
    tryCatch(frm(y ~ x, d, family = family, control = ctl),
             error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, nan = sum(grepl("NA/NaN function evaluation", warns,
                                  fixed = TRUE)),
       other = setdiff(unique(warns), "NA/NaN function evaluation"))
}

cat("(a) and (b): NaN warnings and the optimum, old path against new\n")
cat(sprintf("%-10s %-4s %-4s %-6s %-18s %-18s %-10s %-9s %s\n", "design",
            "old", "new", "mapped", "logLik old", "logLik new",
            "identical", "max|dpar|", "eta at optimum (min)"))
for (nm in names(designs)) {
  d <- designs[[nm]]()
  o <- fit_once(d, "inverse.gaussian", "old")
  w <- fit_once(d, "inverse.gaussian", "new")
  if (inherits(o$fit, "error") || inherits(w$fit, "error")) {
    cat(nm, "ERROR:", conditionMessage(if (inherits(o$fit, "error")) o$fit
                                       else w$fit), "\n")
    next
  }
  same <- identical(o$fit$opt$par, w$fit$opt$par) &&
    identical(as.numeric(logLik(o$fit)), as.numeric(logLik(w$fit)))
  eta <- stats::predict(w$fit, type = "link")
  eta <- if (is.list(eta)) eta[[1]] else eta
  cat(sprintf("%-10s %-4d %-4d %-6s %-18.10f %-18.10f %-10s %-9.2e %.4g\n",
              nm, o$nan, w$nan, format(w$fit$opt$nonfinite_trials),
              as.numeric(logLik(o$fit)), as.numeric(logLik(w$fit)), same,
              max(abs(o$fit$opt$par - w$fit$opt$par)),
              min(as.numeric(unlist(eta)))))
}

cat("\n(c) a genuinely NaN objective still warns\n")
set.seed(20260916)
db <- data.frame(y = rnorm(60), x = rnorm(60))
broken <- frmtmb_family(
  "broken", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  # log of a negative number: NaN at every parameter value, the start
  # included, which is what a mistyped density does
  lpdf = function(y, dpars, aterms) {
    log(-exp(dpars[["sigma"]])) + 0 * dpars[["mu"]]
  })
b <- fit_once(db, broken, "new")
cat("broken family: NaN warnings", b$nan, "; outcome:",
    if (inherits(b$fit, "error")) paste("ERROR:", conditionMessage(b$fit))
    else paste("fit, logLik", format(as.numeric(logLik(b$fit)))), "\n")
