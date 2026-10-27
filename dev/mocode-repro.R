## mo() of integers that do not start at 0. Adapted from the 0.75.0
## brmsnames review's dev/review-brmsnames-mo.R (seed 5, n = 300).
##   Rscript dev/mocode-repro.R base|lane [sample]
## base = rellib-r13 (0.75.0), lane = wt-mocode-lib in front of it.
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$k2 <- d$k + 1L
d$k0 <- d$k - 1L
d$kn <- d$k - 3L
zeta <- c(0.6, 0.3, 0.1)
d$y <- 2 + 1.5 * c(0, cumsum(zeta))[d$k] * 3 + rnorm(n)
q(library(frmtmb))
cat("== frmtmb", format(packageVersion("frmtmb")), arm,
    dirname(system.file(package = "frmtmb")), "\n")
fits <- list()
for (v in c("k0", "k", "k2", "kn")) {
  f <- as.formula(paste0("y ~ mo(", v, ")"))
  w <- character(0)
  fit <- tryCatch(withCallingHandlers(frm(f, data = d),
    warning = function(w0) {
      w <<- c(w, conditionMessage(w0))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("mo(%s): range %s ERROR: %s\n", v,
                paste(range(d[[v]]), collapse = ".."),
                conditionMessage(fit)))
    next
  }
  fits[[v]] <- fit
  fe <- fixef(fit)
  cat(sprintf("mo(%s): range %s logLik %.6f df %d; warnings %d\n", v,
              paste(range(d[[v]]), collapse = ".."),
              as.numeric(logLik(fit)), attr(logLik(fit), "df"), length(w)))
  print(signif(fe[, 1:2], 7))
  zn <- grep("^zeta", names(fit$estimates), value = TRUE)
  cat("simplex parameters:", zn, "lengths",
      lengths(fit$estimates[zn]), "\n")
}
nd <- data.frame(k = 1:4, k0 = 0:3, k2 = 2:5, kn = -2:1)
cat("fitted at the four categories (rows), by fit (columns):\n")
fv <- sapply(fits, function(f) fitted(f, newdata = nd)[, 1])
print(signif(fv, 8))
cat("SE of fitted at the four categories:\n")
print(signif(sapply(fits, function(f) fitted(f, newdata = nd)[, 2]), 6))
cat("fitted(mo(k)) on newdata k = 2..4 only:\n")
print(signif(fitted(fits$k, newdata = data.frame(k = 2:4))[, 1], 8))
for (bad in c(0L, 5L)) {
  r <- tryCatch(fitted(fits$k, newdata = data.frame(k = bad)),
                error = function(e) e)
  cat(sprintf("fitted(mo(k)) at k = %d: %s\n", bad,
              if (inherits(r, "error")) paste("ERROR:", conditionMessage(r))
              else paste(signif(r[, 1], 8))))
}

if (identical(args[2L], "sample") &&
    requireNamespace("frmtmb.sample", quietly = TRUE)) {
  q(library(frmtmb.sample))
  Sys.setenv(FRMTMB_STAN_CACHE = file.path(getwd(), "dev/stan-cache"))
  for (v in c("k0", "k")) {
    w <- character(0)
    ds <- withCallingHandlers(
      frm_sample(fits[[v]], chains = 2, iter = 1000, refresh = 0,
                 seed = 11),
      warning = function(w0) {
        w <<- c(w, conditionMessage(w0))
        invokeRestart("muffleWarning")
      })
    s <- fixef(ds)
    cat(sprintf("frm_sample mo(%s): warnings %d\n", v, length(w)))
    if (length(w)) cat("  ", substr(w, 1, 160), sep = "\n  ")
    print(round(s[, 1:4], 3))
    sv <- as.matrix(ds)[, grep("^simo", variables(ds)), drop = FALSE]
    print(round(colMeans(sv), 3))
  }
}
