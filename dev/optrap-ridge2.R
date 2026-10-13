# Lane optrap: a flat-ridge fixture for test-grad-verdict.R on a response
# with unit 1 (predictor_units() empty), so the base and lane builds take
# the same path. Searches seeds for a fit that converges (code 0) on the
# ridge with the gradient just over grad_tol and a large headroom.
#   Rscript dev/optrap-ridge2.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:200
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
for (s in seeds) {
  set.seed(s)
  n <- 1000
  x1 <- rnorm(n)
  dd <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, 1e-5))
  dd$y <- rnorm(n, 0.5 + 0.5 * dd$x1, 0.8)
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x1 + x2), family = gaussian(),
                               data = dd,
                               control = frmtmb_control(restarts = 0)),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")})
  d <- diagnose(f, quiet = TRUE)
  r1 <- d$max_grad / 1e-3
  r2 <- d$grad_headroom / 1e-3
  if (isTRUE(r1 > 1 && r1 < 100 && r2 > 100)) {
    cat(sprintf("seed %d: unit %d code %d max_grad/tol %.3g headroom/tol %.3g | %s\n",
                s, length(frmtmb:::predictor_units(f$frame)),
                f$opt$convergence, r1, r2,
                paste(substr(w, 1, 70), collapse = " ; ")))
  }
}
