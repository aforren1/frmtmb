# Reviewer of lane optstart: the condition-number branch of the rule.
#   Rscript dev/optstart-rev-kappa.R <lib> <label> [Rhome label]
# 1. kappa of the designs of test-optrap.R's two edited tests, and
#    whether the rule applies there (par_template() against 0.70.0's
#    formula).
# 2. A family of designs x2 = x1 + delta * e, bisected on delta so that
#    the rule's kappa is 1e5 to a few ulps; the branch and the fit from
#    each side, under whichever BLAS this R uses.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), " R", R.home(), "\n")
cat("BLAS:", extSoftVersion()[["BLAS"]], " La:", La_library(), "\n")
ns <- asNamespace("frmtmb")
rule_kappa <- function(fo, d, fam = gaussian()) {
  fr <- frm(fo, data = d, family = fam, dry_run = "frame")
  lp <- fr$linpreds[[1]]
  inp <- ns$predictor_unit_input(fr, lp)
  fit <- ns$predictor_ls(lp, inp$z, inp$ok, inp$w)
  c(rule = fit$kappa, polish = ns$predictor_design_kappa(fr))
}
applies <- function(fo, d) {
  st <- par_template(fo, data = d, family = gaussian())
  !identical(unname(st$betad), log(stats::sd(d$y)))
}
set.seed(6001)
d <- data.frame(x = runif(200, 1, 2))
d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
for (j in 1:6) d[[paste0("x", j)]] <- d$x^j
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6
cat("test-optrap 6001 degree 6: kappa", format(rule_kappa(fo, d)),
    " rule applies:", applies(fo, d), "\n")
set.seed(7017)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7
cat("test-optrap 7017 degree 7: kappa", format(rule_kappa(fo, d)),
    " rule applies:", applies(fo, d), "\n")
# 2. bisection to kappa = 1e5 on the rule's measure
set.seed(42)
n <- 300
base <- data.frame(x1 = rnorm(n), e = rnorm(n), g = factor(rep(1:15, 20)))
base$y <- 2 + base$x1 + rnorm(15, 0, 1)[base$g] + rnorm(n)
mk <- function(delta) {
  d <- base
  d$x2 <- d$x1 + delta * d$e
  d
}
fo2 <- y ~ x1 + x2 + (1 | g)
kap <- function(delta) rule_kappa(fo2, mk(delta))[["rule"]]
lo <- 1e-6; hi <- 1e-3        # kappa(lo) > 1e5 > kappa(hi)
for (i in 1:200) {
  mid <- (lo + hi) / 2
  if (kap(mid) > 1e5) lo <- mid else hi <- mid
  if (hi - lo <= 2 * .Machine$double.eps * hi) break
}
cat(sprintf("delta lo %.17g kappa %.17g applies %s\n", lo, kap(lo),
            applies(fo2, mk(lo))))
cat(sprintf("delta hi %.17g kappa %.17g applies %s\n", hi, kap(hi),
            applies(fo2, mk(hi))))
cat(sprintf("polish kappa at lo %.17g, at hi %.17g\n",
            rule_kappa(fo2, mk(lo))[["polish"]],
            rule_kappa(fo2, mk(hi))[["polish"]]))
for (dl in c(lo, hi)) {
  f <- suppressWarnings(frm(fo2, data = mk(dl)))
  cat(sprintf("fit at delta %.17g: logLik %.12f code %d\n", dl,
              as.numeric(logLik(f)), f$opt$convergence))
}
# the fixed deltas a second R can be asked about
for (dl in c(lo, hi, 0.9999 * lo, 1.0001 * hi)) {
  cat(sprintf("fixed delta %.17g kappa %.17g\n", dl, kap(dl)))
}
if (length(a) > 2) {
  for (dl in as.numeric(strsplit(a[3], ",")[[1]])) {
    cat(sprintf("given delta %.17g kappa %.17g applies %s\n", dl, kap(dl),
                applies(fo2, mk(dl))))
  }
}
