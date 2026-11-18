# rel084-kt copy of dev/reviews/alapse-short-search.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Short review of lane alapse, fix round 1, check 2: the estimated
# constant quantile. For design B (100 rows, y = 1 + 0.5 x + N(0, 1),
# y ~ x; seed 2000 + s) and variants with case weights U(0.2, 3), tied
# data (y rounded to integers, x to one decimal) and both, the global
# maximum of the asymmetric Laplace log-likelihood by brute force:
# S(p) from quantreg::rq.wfit(method = "br") on a grid of 1999 quantiles
# (step 5e-4), the profile W (log p + log(1 - p) - log(S / W) - 1); and
# exactly, from the weighted quantile process (rq.fit.br(tau = -1) on
# the weight-scaled rows), maximized on each interval. Compared with
# the lane fit's own log-likelihood, computed here from its estimates.
#   Rscript dev/reviews/alapse-short-search.R <seeds>
a <- commandArgs(TRUE)
seeds <- eval(parse(text = a[1]))
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(quantreg))
options(frmtmb.notices = FALSE, width = 150)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")

wcls <- function(expr) {
  wc <- character()
  v <- withCallingHandlers(expr, warning = function(z) {
    wc <<- c(wc, class(z)[1]); invokeRestart("muffleWarning")
  }, message = function(z) invokeRestart("muffleMessage"))
  list(v = v, wc = wc)
}
chk <- function(r, p, w) sum(w * r * (p - (r < 0)))
prof <- function(p, S, W) W * (log(p) + log1p(-p) - log(S / W) - 1)

grid_max <- function(X, y, w, grid) {
  W <- sum(w)
  v <- vapply(grid, function(p) {
    f <- suppressWarnings(rq.wfit(X, y, tau = p, weights = w, method = "br"))
    prof(p, chk(y - drop(X %*% f$coefficients), p, w), W)
  }, 0)
  list(ll = max(v), p = grid[which.max(v)])
}
exact_max <- function(X, y, w) {
  W <- sum(w)
  pr <- suppressWarnings(rq.fit.br(X * w, y * w, tau = -1))
  taus <- pr$sol[1, ]
  B <- pr$sol[-(1:3), , drop = FALSE]
  best <- list(ll = -Inf)
  for (j in seq_len(ncol(B))) {
    lo <- max(taus[j], 1e-9)
    hi <- min(if (j < ncol(B)) taus[j + 1] else 1, 1 - 1e-9)
    if (hi <= lo) next
    r <- y - drop(X %*% B[, j])
    pf <- function(p) prof(p, chk(r, p, w), W)
    o <- optimize(pf, c(lo, hi), maximum = TRUE, tol = 1e-12)
    for (pp in c(o$maximum, lo, hi)) {
      v <- pf(pp)
      if (is.finite(v) && v > best$ll) best <- list(ll = v, p = pp)
    }
  }
  best
}

grid <- seq(5e-4, 1 - 5e-4, by = 5e-4)
res <- NULL
for (s in seeds) for (vt in c("plain", "w", "ties", "wties")) {
  set.seed(2000 + s)
  d <- data.frame(x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  d$w <- 1
  if (vt %in% c("w", "wties")) d$w <- runif(100, 0.2, 3)
  if (vt %in% c("ties", "wties")) {
    d$y <- round(d$y)
    d$x <- round(d$x, 1)
  }
  r <- if (vt %in% c("w", "wties")) {
    wcls(frm(y | weights(w) ~ x, family = asym_laplace(), data = d))
  } else wcls(frm(y ~ x, family = asym_laplace(), data = d))
  f <- r$v
  X <- cbind(i = 1, x = d$x)
  dp <- ns$eval_dpars(f)[[1]]
  b <- ns$fixef_by_dpar(f)$mu
  p <- dp$quantile[1]
  sg <- dp$sigma[1]
  rr <- d$y - drop(X %*% b)
  ll <- sum(d$w * (log(p) + log1p(-p) - log(sg) - rr * (p - (rr < 0)) / sg))
  gm <- grid_max(X, d$y, d$w, grid)
  ex <- exact_max(X, d$y, d$w)
  cert <- f$opt$ald_cert
  cat(sprintf(paste0("%-5s %3d code %d global %-5s warn %-40s ll %.6f ",
                     "grid %.6f exact %.6f gap_grid %+.2e gap_exact %+.2e ",
                     "p %.4f p_grid %.4f p_exact %.4f evals %s moved %s\n"),
              vt, s, f$opt$convergence, isTRUE(cert$global),
              paste(r$wc, collapse = ","), ll, gm$ll, ex$ll, gm$ll - ll,
              ex$ll - ll, p, gm$p, ex$p, cert$search$evals %||% NA,
              format(cert$search$moved_from %||% NA, digits = 4)))
  res <- rbind(res, data.frame(vt = vt, s = s, code = f$opt$convergence,
                               global = isTRUE(cert$global),
                               ald_local = "frmtmb_ald_local" %in% r$wc,
                               conv = any(r$wc %in% c("frmtmb_conv_code",
                                                      "frmtmb_conv_gradient")),
                               gg = gm$ll - ll, ge = ex$ll - ll))
}
for (vt in unique(res$vt)) {
  z <- res[res$vt == vt, ]
  cat(sprintf(paste0("SUMMARY %-5s n %d | code0 %d | global %d | ald_local %d",
                     " | conv-warn %d | gap_grid > 1e-6: %d (max %.3g) | ",
                     "gap_exact > 1e-6: %d (max %.3g) | global & gap_exact > 1e-6: %d\n"),
              vt, nrow(z), sum(z$code == 0), sum(z$global), sum(z$ald_local),
              sum(z$conv), sum(z$gg > 1e-6), max(z$gg), sum(z$ge > 1e-6),
              max(z$ge), sum(z$global & z$ge > 1e-6)))
}
cat("DONE\n")
