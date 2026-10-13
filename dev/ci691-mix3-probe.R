# Lane ci691: the three-component ordinal mixture of
# test-nonfinite-gradient.R ("the three-component mixture the review
# found ends finite"), whose fit let "singular convergence (7)" escape
# on the ubuntu runners. Fits it as the test does (k = 0) and with each
# x moved by a few ulps (k > 0), a stand-in for another platform's
# rounding, and prints the code, message, objective, whether the
# gradient is finite, and the warnings.
#   Rscript dev/ci691-mix3-probe.R <lib> <kmax>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
kmax <- if (length(args) > 1) as.integer(args[2]) else 12L
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(20261005 + 31)
n <- 400
x <- rnorm(n)
invisible(rnorm(n))
invisible(sample(c("a", "b"), n, TRUE))
invisible(rbinom(n, 1, 0.4))
invisible(rlogis(n))
invisible(runif(n))
invisible(runif(n, 0.5, 2))
cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
d <- data.frame(y = as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf))),
                x = x)
fam <- mixture(cumulative("probit"), sratio("cloglog"), acat())
u <- suppressWarnings(frm(bf(y ~ x), family = fam, data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
for (k in 0:kmax) {
  dk <- d
  if (k > 0) {
    set.seed(k)
    dk$x <- d$x * (1 + sample(-4:4, n, TRUE) * .Machine$double.eps)
  }
  ws <- character(0)
  fit <- withCallingHandlers(
    frm(bf(y ~ x), family = fam, data = dk),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  o <- fit$opt
  cat(sprintf("k %2d code %d %-28s obj %.10f grad_finite %s\n", k,
              o$convergence, o$message, o$objective,
              all(is.finite(fit$obj$gr(o$par)))))
  for (w in ws) cat("     ", substr(w, 1, 90), "\n")
}
