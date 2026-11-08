# Review lane osa: randomization, RNG state, and an independent poisson
# reference built from ppois()/dpois() at mu = exp(b0 + b1 x) from the
# printed coefficients (no frmtmb internals).
#
#   Rscript dev/review-osa-rng.R [lib]
args <- commandArgs(TRUE)
lib1 <- if (length(args) >= 1) args[1] else "C:/Users/adf44/source/r/wt-osa-lib"
.libPaths(c(if (lib1 != "none") lib1, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-rng\n")
# data built without touching the global RNG stream: a local generator
d <- local({
  x <- qnorm(ppoints(150))
  x <- x[order(sin(seq_along(x) * 7.1))]
  mu <- exp(0.8 + 0.5 * x)
  # deterministic counts from the quantile function, no RNG
  y <- qpois(ppoints(150)[order(cos(seq_along(x) * 3.3))], mu)
  data.frame(x = x, y = y)
})
fit <- frm(y ~ x, family = poisson(), data = d)
cat("1. .Random.seed exists after fit, before residuals:",
    exists(".Random.seed", envir = globalenv(), inherits = FALSE), "\n")
r1 <- residuals(fit, type = "osa")[, "Estimate"]
ex <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
cat("2. after residuals on a session with no seed: exists", ex,
    " class", if (ex) class(get(".Random.seed", envir = globalenv())) else
      "-", "\n")
w <- NULL
withCallingHandlers(invisible(runif(1)), warning = function(c) {
  w <<- conditionMessage(c)
  invokeRestart("muffleWarning")
})
cat("3. next runif(1) warning:", if (is.null(w)) "none" else w, "\n")

set.seed(42)
s0 <- .Random.seed
r2 <- residuals(fit, type = "osa")[, "Estimate"]
cat("4. seeded session: .Random.seed unchanged by residuals():",
    identical(s0, .Random.seed), "\n")
cat("5. two calls identical:", identical(r1, r2), "\n")
r3 <- residuals(fit, type = "osa", seed = 7)[, "Estimate"]
cat("6. seed = 7 passed through differs from default:",
    !identical(r1, r3), " max|diff|", format(max(abs(r1 - r3)), digits = 3),
    "\n")
set.seed(42)
s0 <- .Random.seed
r4 <- tryCatch(residuals(fit, type = "osa", seed = NULL)[, "Estimate"],
               error = function(e) e)
if (inherits(r4, "error")) {
  cat("7. seed = NULL:", conditionMessage(r4), "\n")
} else {
  cat("7. seed = NULL advances the session RNG:", !identical(s0, .Random.seed),
      "\n")
}

# independent reference: TMB's documented default uniforms
U <- local({ set.seed(123); runif(nrow(d)) })
b <- fixef(fit)[, "Estimate"]
mu <- exp(b[[1]] + b[[2]] * d$x)
ref <- qnorm(ppois(d$y, mu) - U * dpois(d$y, mu))
cat("8. poisson vs qnorm(ppois(y) - U dpois(y)), U = set.seed(123):",
    "max|diff|", format(max(abs(r1 - ref)), digits = 3), "\n")
# same but against the other natural orientation, to show the test can fail
ref2 <- qnorm(ppois(d$y - 1, mu) + U * dpois(d$y, mu))
cat("   against qnorm(F(y-1) + U p(y)) instead:",
    format(max(abs(r1 - ref2)), digits = 3), "\n")
# a dataset with a different response gets the same U sequence
d2 <- d
d2$y <- rev(d$y)
f2 <- frm(y ~ x, family = poisson(), data = d2)
b2 <- fixef(f2)[, "Estimate"]
mu2 <- exp(b2[[1]] + b2[[2]] * d2$x)
r5 <- residuals(f2, type = "osa")[, "Estimate"]
u5 <- (ppois(d2$y, mu2) - pnorm(r5)) / dpois(d2$y, mu2)
cat("9. second data set reads the same uniforms: max|u - U|",
    format(max(abs(u5 - U)), digits = 3), "\n")
# subset of rows: is U indexed by row or by position?
r6 <- residuals(fit, type = "osa", subset = 76:150)
cat("10. subset = 76:150 rows returned:", NROW(r6), "\n")
cat("END review-osa-rng\n")
