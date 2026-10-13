# Lane optrap: test-se-tier.R's RQ1 case (quadrature = TRUE on a fit
# left short) rebuilt by its start, on a library.
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(36)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
d$y <- 1 + 0.5 * d$x + rnorm(100)
l4 <- suppressMessages(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
for (st in list(NULL, list(theta = -20))) {
  w <- character(); m <- list()
  f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d, quadrature = TRUE,
                               start = st),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m[[length(m) + 1L]] <<- x
      invokeRestart("muffleMessage")})
  cat(sprintf("start %s: logLik - lme4 %.4f | boundary messages %d | se_lost %s\n",
              if (is.null(st)) "default" else "theta = -20",
              as.numeric(logLik(f)) - as.numeric(logLik(l4)),
              sum(vapply(m, inherits, NA, "frmtmb_boundary_fit")),
              paste(frmtmb:::sdr_of(f)$se_lost, collapse = ",")))
  for (x in w) cat("   W:", substr(x, 1, 120), "\n")
}
