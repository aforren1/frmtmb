# The one fixture whose optimizer evaluated a floored log sd below
# -308.45 (dev/moparity-audit-run3.R: test-predfix.R, the (1 + x | g)
# block of "a random slope on a tiny column reaches the scale-1
# optimum"): its fits with the floor and, in the same process, with
# sd_floored() bound to exp().
#   Rscript dev/moparity-predfix-floor.R
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
reg_new <- get("covstruct_registry", ns)
reg_old <- reg_new
old_env <- new.env(parent = ns)
old_env$sd_floored <- exp
for (nm in names(reg_old)) {
  f <- reg_old[[nm]]$nll
  if (is.function(f) && identical(environment(f), ns)) {
    environment(f) <- old_env
    reg_old[[nm]]$nll <- f
  }
}
use <- function(reg) {
  unlockBinding("covstruct_registry", ns)
  assign("covstruct_registry", reg, envir = ns)
  lockBinding("covstruct_registry", ns)
}
set.seed(2)
g <- factor(rep(1:20, each = 15))
x <- stats::rnorm(300)
d <- data.frame(g = g, x0 = x,
                y = 1 + 0.5 * x + stats::rnorm(20, 0, 0.7)[g] +
                  stats::rnorm(20, 0, 0.4)[g] * x + stats::rnorm(300))
fo <- bf(y ~ x + (1 + x | g))
for (s in c(1, 1e-2, 1e-3, 1e-6)) {
  d$x <- d$x0 * s
  res <- list()
  for (arm in c("floor", "exp")) {
    use(if (arm == "floor") reg_new else reg_old)
    f <- suppressWarnings(frm(fo, family = gaussian(), data = d))
    res[[arm]] <- list(ll = as.numeric(logLik(f)), par = f$opt$par,
                       evals = f$opt$evals)
  }
  cat(sprintf("scale %g: logLik floor %.12f exp %.12f | estimates identical %s | logLik identical %s\n",
              s, res$floor$ll, res$exp$ll,
              identical(res$floor$par, res$exp$par),
              identical(res$floor$ll, res$exp$ll)))
}
# the same file's "a pre-fit that errors falls back under the default"
set.seed(533)
g <- factor(rep(1:20, each = 12))
d <- data.frame(g, x = stats::rnorm(240) * 0.03, z = stats::rnorm(240))
d$yb <- as.integer(d$z > 0)
fo <- yb ~ z + x + (1 + x | g)
for (ctl in list(list(), list(autoscale = FALSE), list(autoscale = TRUE))) {
  res <- list()
  for (arm in c("floor", "exp")) {
    use(if (arm == "floor") reg_new else reg_old)
    f <- tryCatch(suppressWarnings(frm(fo, family = bernoulli(), data = d,
                                       control = do.call(frmtmb_control,
                                                         ctl))),
                  error = function(e) conditionMessage(e))
    res[[arm]] <- if (is.character(f)) f else f$opt$par
  }
  cat(sprintf("bernoulli seed 533, control %s: identical %s (%s)\n",
              paste(names(ctl), unlist(ctl), collapse = " "),
              identical(res$floor, res$exp),
              if (is.character(res$floor)) paste("error:",
                substr(res$floor, 1, 60)) else "a fit"))
}
# what each arm says with autoscale = TRUE, in full
for (arm in c("floor", "exp")) {
  use(if (arm == "floor") reg_new else reg_old)
  w <- character(0)
  f <- tryCatch(withCallingHandlers(
    frm(fo, family = bernoulli(), data = d,
        control = frmtmb_control(autoscale = TRUE)),
    warning = function(z) {
      w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(f, "error")) {
    cat(arm, ": error:", conditionMessage(f), "\n")
  } else {
    cat(arm, ": a fit, logLik", format(as.numeric(logLik(f)), digits = 12),
        "code", f$opt$convergence, "| sd", format(exp(f$estimates$theta[1:2]),
                                                  digits = 4),
        "| warnings:", paste(substr(w, 1, 120), collapse = " || "), "\n")
  }
}
