# rel070 seflags: attribution of the false-alarm battery's differences.
# Rebuilds the battery's poly and bounded fits (dev/seflags-rev3-falarm.R,
# same data and the same order of fits) and prints, for the fits named,
# the optimizer code, logLik, every warning, the lost set, the outer SEs,
# predictor_units() (optrap's response units) and opt$polish (optrap's
# Newton polish) where the build has them.
#   Rscript rel070-seflags-attrib.R <lib> <seeds expr>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "| R", R.home(), "\n")
seeds <- eval(parse(text = a[2]))
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
show <- function(tag, r) {
  f <- r$v
  l <- suppressWarnings(ns$sdr_of(f)$se_lost)
  u <- if (exists("predictor_units", ns)) ns$predictor_units(f$frame) else NA
  se <- suppressWarnings(sqrt(diag(as.matrix(vcov(f, full = TRUE)))))
  cat(sprintf("%s | code %s | logLik %.10g | lost %s | units %s | polish %s\n",
              tag, f$opt$convergence, as.numeric(logLik(f)),
              if (length(l)) paste(names(l), l, sep = ":", collapse = " ") else "-",
              if (length(u)) paste(names(u), signif(u, 4), collapse = " ") else "-",
              if (is.null(f$opt$polish)) "-" else
                paste(names(f$opt$polish), signif(f$opt$polish, 4),
                      collapse = " ")))
  cat("   par", paste(signif(f$opt$par, 8), collapse = " "), "\n")
  cat("   se ", paste(signif(se, 5), collapse = " "), "\n")
  for (w in r$w) cat("   W:", substr(gsub("\n", " ", w), 1, 300), "\n")
}
for (s in seeds) {
  set.seed(s)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
  for (deg in 5:6) {
    rhs <- paste(paste0("x", 1:deg), collapse = " + ")
    show(sprintf("seed %d poly%d_lm  ", s, deg),
         cap(frm(stats::as.formula(paste("y ~", rhs)), data = d)))
    show(sprintf("seed %d poly%d_lmer", s, deg),
         cap(frm(stats::as.formula(paste("y ~", rhs, "+ (1 | g)")),
                 data = d)))
  }
  set.seed(s)
  n <- 150
  x <- rnorm(n)
  yy <- 0.2 + 0.5 * x + rnorm(n)
  mu <- plogis(-2 + 0.5 * x)
  yb <- rbeta(n, mu * 200, (1 - mu) * 200)
  yb[sample(n, 3)] <- 0
  yb2 <- rbeta(n, 0.98 * 20, 0.02 * 20)
  yb2[sample(n, 2)] <- 1
  yt <- 1 + 0.5 * x + rt(n, df = 30)
  dn <- data.frame(x = x, c = rpois(n, exp(1 + 0.3 * x)))
  r <- cap(frm(c ~ x, family = negbinomial(), data = dn))
  show(sprintf("seed %d negbin_large_shape", s), r)
  he <- r$v$obj$he(r$v$opt$par)
  cat("   exact-Hessian SE", paste(signif(suppressWarnings(
    sqrt(diag(solve(he)))), 5), collapse = " "), "| eigen",
    paste(signif(eigen(he, only.values = TRUE)$values, 4), collapse = " "),
    "\n")
}
