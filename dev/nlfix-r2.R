# R2: no stop on a correct model. Exact zero slope, slope 2e-5, null data sets, and
# the sharp maximum of nearly noiseless data, on the review's seeds.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
tt <- function(expr) {
  msg <- character(0); warn <- list()
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warn[[length(warn) + 1L]] <<- w; invokeRestart("muffleWarning") })
  list(r = r, msg = msg, warn = warn)
}
named <- function(a) Filter(function(w) inherits(w, "frmtmb_nl_saddle"), a$warn)
fb <- bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE)
set.seed(11003)
n <- 300
dz <- data.frame(x = rnorm(n)); dz$x <- dz$x - mean(dz$x)
e <- rnorm(n)
dz$y0 <- 1 + as.numeric(resid(lm(e ~ x, dz)))
for (s in c(0, 2e-5)) {
  dz$y <- dz$y0 + s * dz$x
  a <- tt(frm(fb, data = dz))
  cat(sprintf("slope %g: %s ll %.5f, linear %.5f, code %s, saddle warnings %d\n", s,
              if (inherits(a$r, "error")) "ERROR" else "fit",
              as.numeric(logLik(a$r)), as.numeric(logLik(frm(y ~ x, data = dz))),
              a$r$opt$convergence, length(named(a))))
  for (w in named(a)) cat("  WARN:", conditionMessage(w), "\n")
}
set.seed(11002)
nerr <- 0L; nw <- 0L
for (i in 1:400) {
  d0 <- data.frame(x = rnorm(40), y = rnorm(40))
  a <- tt(frm(fb, data = d0))
  if (inherits(a$r, "error")) nerr <- nerr + 1L else if (length(named(a))) nw <- nw + 1L
}
cat("null data sets, 400: errors", nerr, " saddle warnings", nw, "\n")
for (sg in c(1e-3, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)) {
  set.seed(11004)
  dn <- data.frame(x = rnorm(200)); dn$y <- 1 + 2 * dn$x + rnorm(200, 0, sg)
  lin <- as.numeric(logLik(frm(y ~ x, data = dn)))
  a <- tt(frm(fb, data = dn))
  b <- tt(frm(fb, data = dn, start = list(beta = c("a_(Intercept)" = 1, "b1_(Intercept)" = 1, "b2_(Intercept)" = 1))))
  cat(sprintf("sigma %.0e: linear %.4f | default %s code %s, saddle warnings %d | started %s code %s\n", sg, lin,
              if (inherits(a$r, "error")) "ERROR" else format(as.numeric(logLik(a$r)), digits = 9),
              if (inherits(a$r, "error")) "NA" else a$r$opt$convergence, length(named(a)),
              if (inherits(b$r, "error")) "ERROR" else format(as.numeric(logLik(b$r)), digits = 9),
              if (inherits(b$r, "error")) "NA" else b$r$opt$convergence))
}
