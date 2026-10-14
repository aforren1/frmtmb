# Review of lane ci701: the degenerate ordinal mixture of
# test-ordinal-mixture.R:359 on ulp-moved x (seed 9000 + k, not the
# lane's 2000 + k). For each fit: the code, the log-likelihood, and how
# much an independent polish gains: optim(BFGS) from the estimate, and
# nlminb restarted from the estimate moved by 0.05 in every coordinate.
#
#   Rscript dev/ci701-rev-omx.R <k0> <k1>
args <- commandArgs(TRUE)
k0 <- as.integer(args[1]); k1 <- as.integer(args[2])
.libPaths(c("C:/Users/adf44/source/r/wt-ci701-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < stats::plogis(-1.2 + 0.5 * d$z), 0L, d$y)
  d
}
d0 <- omx_data(20261011, n = 300)
for (k in k0:k1) {
  d <- d0
  if (k > 0) {
    set.seed(9000 + k)
    d$x <- d$x * (1 + sample(-3:3, nrow(d), TRUE) * 2^-52)
  }
  w <- character(0)
  fit <- withCallingHandlers(
    frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(), sratio()),
        data = d),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  obj <- fit$obj
  p <- fit$opt$par
  f0 <- obj$fn(p)
  ok <- function(v) if (is.finite(v)) v else Inf
  o1 <- tryCatch(stats::optim(p, function(q) ok(obj$fn(q)), obj$gr,
                              method = "BFGS",
                              control = list(maxit = 2000, reltol = 1e-14)),
                 error = function(e) list(value = NA))
  set.seed(k)
  o2 <- tryCatch(stats::nlminb(p + 0.05 * sign(stats::rnorm(length(p))),
                               function(q) ok(obj$fn(q)), obj$gr,
                               control = list(eval.max = 3000,
                                              iter.max = 2000)),
                 error = function(e) list(objective = NA))
  cat(sprintf(paste0("RES k=%d code=%d msg=%s ll=%.8f gain_bfgs=%.3g ",
                     "gain_restart=%.3g nwarn=%d w=%s\n"), k,
              fit$opt$convergence, gsub(" ", "_", fit$opt$message), -f0,
              f0 - o1$value, f0 - o2$objective, length(w),
              paste(substr(w, 1, 40), collapse = " | ")))
}
