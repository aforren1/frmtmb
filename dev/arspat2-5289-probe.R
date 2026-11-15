# Backlog 5289: student() with fcor(), mi() and something else integrated
# (a random effect, beta under REML or profile). The ML objective of the
# mi() fit against the same model fit on the observed rows only
# (fcor(M_oo), the margin), at the same outer parameters: one model, so
# the two objectives must agree. Base build when ARSPAT2_LIB="".
# Data seed 4 (as dev/reviews/arspat-short-a2.R), n = 60, m = 20.
lib <- Sys.getenv("ARSPAT2_LIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
ns <- asNamespace("frmtmb")
n <- 60; m <- 20
M <- 0.5^abs(outer(1:n, 1:n, "-"))
mk <- function(re_sd, seed = 4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)))
  d$y <- 1 + 0.5 * d$x + rnorm(6, 0, re_sd)[d$h] +
    as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
  miss <- sort(sample(n, m)); d$y[miss] <- NA
  d
}
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
cases <- list(
  list(lab = "(1 | h), re_sd 0", f = y ~ x + (1 | h) + fcor(M), sd = 0),
  list(lab = "(1 | h), re_sd 1.5", f = y ~ x + (1 | h) + fcor(M), sd = 1.5),
  list(lab = "no RE (constant route)", f = y ~ x + fcor(M), sd = 0),
  list(lab = "no RE, REML", f = y ~ x + fcor(M), sd = 0, REML = TRUE),
  list(lab = "(1 | h), profile", f = y ~ x + (1 | h) + fcor(M), sd = 1.5,
       profile = TRUE)
)
worst <- 0
for (cs in cases) {
  d <- mk(cs$sd)
  ok <- !is.na(d$y)
  fm <- update(cs$f, y | mi() ~ .)
  fo <- cs$f
  fo[[3]] <- do.call(substitute, list(fo[[3]], list(M = quote(Mo))))
  ctl <- frmtmb_control(profile = isTRUE(cs$profile))
  reml <- isTRUE(cs$REML)
  f <- quiet(frm(fm, data = d, data2 = list(M = M), family = student(),
                 REML = reml, control = ctl))
  g <- quiet(frm(fo, data = d[ok, ], data2 = list(Mo = M[ok, ok]),
                 family = student(), REML = reml, control = ctl))
  lab <- ns$outer_par_labels(f)
  cat(sprintf("\n=== %s: labels %s\n", cs$lab, paste(lab, collapse = " ")))
  inu <- grep("^nu", lab)
  # the same outer point in both objectives, at nu near 3, 6 and 21
  for (k in 0:2) {
    p <- f$opt$par + 0.05 * k * seq_along(f$opt$par) / length(f$opt$par)
    p[inu] <- c(log(2), log(5), log(20))[k + 1L]
    a <- -f$obj$fn(p)
    b <- -g$obj$fn(p)
    worst <- max(worst, abs(a - b))
    cat(sprintf("  point %d: -fn mi %.8f  -fn observed rows %.8f  diff %+.3e\n",
                k, a, b, a - b))
  }
  cat("  mi_laplace_exact", ns$mi_laplace_exact(f$frame), " tk",
      format(f$mi_laplace$tk %||% NA), "\n")
}
cat(sprintf("\nlargest |diff| over every point: %.3e\n", worst))
