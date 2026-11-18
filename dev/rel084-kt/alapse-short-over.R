# rel084-kt copy of dev/reviews/alapse-short-over.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Short review of lane alapse: frmtmb_ald_local on fits whose maximum is
# provably the highest by the lane's own argument (a constant sigma on
# another link is the same reparameterization; two responses that share
# no parameter separate). Compared with the proved-global fits.
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cap <- function(expr) {
  wc <- character()
  v <- withCallingHandlers(expr, warning = function(z) {
    wc <<- c(wc, class(z)[1]); invokeRestart("muffleWarning")
  }, message = function(z) invokeRestart("muffleMessage"))
  list(v = v, wc = wc)
}
set.seed(7)
d <- data.frame(x = rnorm(200), z = rnorm(200))
d$y <- 1 + 0.5 * d$x + exp(0.4 * d$z) * rnorm(200)
d$y2 <- -1 + 0.3 * d$x + rnorm(200)
ll <- function(r) sprintf("logLik %.8f code %d warn %s",
                          as.numeric(logLik(r$v)), r$v$opt$convergence,
                          paste(r$wc, collapse = ","))
a <- cap(frm(bf(y ~ x, quantile = 0.3), family = asym_laplace(), data = d))
b <- cap(frm(bf(y ~ x, quantile = 0.3),
             family = asym_laplace(link_sigma = "softplus"), data = d))
cat("log link     :", ll(a), "\n")
cat("softplus link:", ll(b), "\n")
cat(sprintf("diff %.3g\n", as.numeric(logLik(a$v)) - as.numeric(logLik(b$v))))
c2 <- cap(frm(bf(y2 ~ x, quantile = 0.5), family = asym_laplace(), data = d))
a5 <- cap(frm(bf(y ~ x, quantile = 0.5), family = asym_laplace(), data = d))
m <- cap(frm(bf(y ~ x, quantile = 0.5) + bf(y2 ~ x, quantile = 0.5),
             family = asym_laplace(), data = d))
cat("y alone      :", ll(a5), "\n")
cat("y2 alone     :", ll(c2), "\n")
cat("both         :", ll(m), "\n")
cat(sprintf("sum of singles - joint %.3g\n",
            as.numeric(logLik(a5$v)) + as.numeric(logLik(c2$v)) -
              as.numeric(logLik(m$v))))
