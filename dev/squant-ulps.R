# Lane squant: do the new tests' fit assertions hold on every rounding?
# The response of each test design times (1 + 2e-15 k), k = 0..20, and
# the assertion re-evaluated: (a) the random-intercept fit stops at
# code 0 with no warning and both profile points at log sd -2 and +2
# are below it; (b) the gap to rq() shrinks with the smoothing and is
# below 0.25 standard errors at 0.02; (c) the estimating equation holds
# to 1e-4.
#
# Run: Rscript dev/squant-ulps.R > dev/squant-out/ulps.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "BLAS", extSoftVersion()[["BLAS"]], "\n")
ns <- asNamespace("frmtmb")
nwarn <- function(expr) {
  ws <- 0L
  v <- withCallingHandlers(expr, warning = function(w) {
    ws <<- ws + 1L
    invokeRestart("muffleWarning")
  })
  list(v = v, w = ws)
}

set.seed(20261114)
n <- 80
dre <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 4), each = 4)))
dre$y <- 1 + 0.5 * dre$x + rnorm(n / 4, 0, 0.7)[as.integer(dre$g)] + rnorm(n)
set.seed(41)
drq <- data.frame(x = runif(300))
drq$y <- 1 + 2 * drq$x + (0.5 + drq$x) * rnorm(300)
set.seed(42)
dee <- data.frame(x = rnorm(400))
dee$y <- 0.5 + dee$x + rexp(400)
r <- quantreg::rq(y ~ x, tau = 0.75, data = drq)
se <- summary(r, se = "nid")$coefficients[, 2]

ok <- c(a = 0L, b = 0L, c = 0L)
for (k in 0:20) {
  d <- dre
  d$y <- d$y * (1 + 2e-15 * k)
  f <- nwarn(frm(bf(y ~ x + (1 | g), quantile = 0.5),
                 family = quantile_elf(smoothing = 0.5), data = d))
  fit <- f$v
  par <- fit$opt$par
  it <- which(names(par) == "theta")[1L]
  prof <- vapply(c(-2, 2), function(s) {
    -nlminb(par[-it], function(p) {
      full <- par
      full[-it] <- p
      full[it] <- par[it] + s
      fit$obj$fn(full)
    })$objective
  }, 0)
  a <- fit$opt$convergence == 0 && f$w == 0 && all(prof < -fit$opt$objective)
  d <- drq
  d$y <- d$y * (1 + 2e-15 * k)
  gap <- vapply(c(1, 0.25, 0.02), function(lam) {
    ff <- nwarn(frm(bf(y ~ x, quantile = 0.75),
                    family = quantile_elf(smoothing = lam), data = d))
    if (ff$w) return(NA)
    max(abs(fixef(ff$v)[, "Estimate"] - coef(r)) / se)
  }, 0)
  b <- all(is.finite(gap)) && gap[3] < gap[1] && gap[3] < 0.25
  d <- dee
  d$y <- d$y * (1 + 2e-15 * k)
  cc <- TRUE
  for (p in c(0.2, 0.9)) {
    ff <- nwarn(frm(bf(y ~ x, quantile = p), family = quantile_elf(),
                    data = d))
    dp <- lapply(ns$eval_dpars(ff$v)[[1L]], as.numeric)
    t <- (d$y - dp$mu) / (0.25 * dp$sigma)
    e1 <- abs(mean(plogis(t)) - (1 - p)) / (1 - p)
    e2 <- abs(mean(d$x * (plogis(t) - (1 - p)))) / mean(abs(d$x))
    cc <- cc && ff$w == 0 && e1 < 1e-4 && e2 < 1e-4
  }
  ok <- ok + c(a, b, cc)
  cat(sprintf(paste0("k %2d: (a) %s conv %d warn %d profile below by %.3f",
                     " %.3f | (b) %s gaps %s | (c) %s\n"), k, a,
              fit$opt$convergence, f$w, -fit$opt$objective - prof[1],
              -fit$opt$objective - prof[2], b,
              paste(sprintf("%.3f", gap), collapse = " "), cc))
}
cat("HELD of 21: (a)", ok[["a"]], "(b)", ok[["b"]], "(c)", ok[["c"]], "\n")
