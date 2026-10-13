# Lane optrap: the flat-ridge fixture of test-grad-verdict.R at two
# scales of the response, on a library.
#   Rscript dev/optrap-ridge.R <lib or "">
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
for (k in c(1, 0.5, 0.25)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  dd <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, 1e-5))
  dd$y <- rnorm(n, 1 + 2 * dd$x1, 1) * k
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x1 + x2), family = gaussian(),
                               data = dd,
                               control = frmtmb_control(restarts = 0)),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")})
  d <- diagnose(f, quiet = TRUE)
  cat(sprintf("k %g: code %d max_grad/tol %.3g headroom/tol %.3g | %s\n", k,
              f$opt$convergence, d$max_grad / 1e-3, d$grad_headroom / 1e-3,
              paste(substr(w, 1, 60), collapse = " ; ")))
}
# a start on the ridge, far along it
for (dist in c(1e4, 3e4)) for (it in c(1, 1000)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  dd <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, 1e-5))
  dd$y <- rnorm(n, 1 + 2 * dd$x1, 1)
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ x1 + x2), family = gaussian(),
                               data = dd, start = list(beta = c(1, 2 + dist, -dist), betad = 0),
                               control = frmtmb_control(restarts = 0,
                                 optCtrl = list(iter.max = it, eval.max = 10))),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")})
  d <- diagnose(f, quiet = TRUE)
  cat(sprintf("dist %g it %d: code %d max_grad/tol %.3g headroom/tol %.3g | %s\n",
              dist, it, f$opt$convergence, d$max_grad / 1e-3,
              d$grad_headroom / 1e-3, paste(substr(w, 1, 60), collapse = " ; ")))
}
