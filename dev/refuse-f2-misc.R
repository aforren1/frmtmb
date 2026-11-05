# N5 calibration, the review's second table (dev/review-refuse-f1-misc.R, seed
# 11501, the same draws): exgaussian, reader slope r, 20 or 60 read of 200. For
# each cell: the shift of y's dpars from the observed rows, in SE, and what the
# warning says. Then the cost of the check as a share of a fit.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
pr <- function(...) cat(sprintf(...), "\n", sep = "")
fit_w <- function(e) {
  w <- character()
  v <- withCallingHandlers(e, warning = function(x) {
    w <<- c(w, class(x)[1L]); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(fit = v, warned = "frmtmb_mi_laplace_error" %in% w)
}
# frm()'s shift from an exact maximum by quadrature over each missing value,
# as dev/review-refuse-f1-laplace.R writes it: no frmtmb density in it
exact_exg <- function(d, fit) {
  dexg <- function(y, mu, sig, bet) {
    m <- mu - bet
    -log(bet) + (m - y) / bet + sig^2 / (2 * bet^2) +
      pnorm((y - m) / sig - sig / bet, log.p = TRUE)
  }
  rows <- which(is.na(d$y)); ob <- setdiff(seq_len(nrow(d)), rows)
  g <- seq(-15, 40, length.out = 1401); dg <- diff(g)[1]
  nll <- function(p) {
    ey <- p[1] + p[2] * d$x
    s <- sum(dexg(d$y[ob], ey[ob], exp(p[3]), exp(p[4]))) +
      sum(dnorm(d$z[ob], p[5] + p[6] * d$y[ob], exp(p[7]), log = TRUE))
    for (i in rows) {
      li <- dexg(g, ey[i], exp(p[3]), exp(p[4])) +
        dnorm(d$z[i], p[5] + p[6] * g, exp(p[7]), log = TRUE)
      li[!is.finite(li)] <- -Inf
      mx <- max(li); s <- s + mx + log(sum(exp(li - mx)) * dg)
    }
    if (is.finite(s)) -s else 1e10
  }
  e <- fit$estimates
  p0 <- c(e$beta[grep("^y_", names(e$beta))], e$betad[grep("^y_", names(e$betad))],
          e$beta[grep("^z_", names(e$beta))], e$betad[grep("^z_", names(e$betad))])
  o <- nlminb(p0, nll, control = list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-10))
  se <- sqrt(diag(solve(stats::optimHess(o$par, nll))))
  max(abs((p0 - o$par) / se))
}
set.seed(11501)
for (r in c(0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7)) for (mm in c(20, 60)) {
  shs <- tot <- wn <- numeric(8)
  for (k in 1:8) {
    n <- 200; x <- rnorm(n)
    y <- 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2
    d <- data.frame(x = x, y = y); d$z <- 0.5 + r * y + rnorm(n, 0, 1)
    rows <- sample(n, mm); d$y[rows] <- NA
    j <- tryCatch(fit_w(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y)) + gaussian(), data = d)),
                  error = function(e) NULL)
    cc <- tryCatch(suppressWarnings(frm(bf(y ~ x), family = exgaussian(), data = d[-rows, ])),
                   error = function(e) NULL)
    if (is.null(j) || is.null(cc)) { shs[k] <- tot[k] <- wn[k] <- NA; next }
    H <- stats::optimHess(cc$opt$par, cc$obj$fn, cc$obj$gr); se <- sqrt(diag(solve(H)))
    ej <- j$fit$estimates$betad[grep("^y_", names(j$fit$estimates$betad))]
    shs[k] <- max(abs((ej - cc$estimates$betad) / se[3:4]))
    tot[k] <- j$fit$mi_laplace$total
    wn[k] <- j$warned
    # a silent row with a shift past 1 SE: from the observed rows, a reader
    # that carries information is not an error, so read it against the exact
    # maximum, which has no Laplace in it
    if (!j$warned && shs[k] > 1) {
      ex <- exact_exg(d, j$fit)
      pr("SILENT slope %.2f read %2d rep %d: shift from observed %.2f, total %.2f, largest shift from the exact maximum %.2f",
         r, mm, k, shs[k], tot[k], ex)
    }
  }
  pr("MISC slope %.2f read %2d: shift SE median %.2f max %.2f; |total| median %.2f min %.2f; warned %d of %d",
     r, mm, median(shs, na.rm = TRUE), max(shs, na.rm = TRUE),
     median(abs(tot), na.rm = TRUE), min(abs(tot), na.rm = TRUE),
     sum(wn, na.rm = TRUE), sum(!is.na(wn)))
}

