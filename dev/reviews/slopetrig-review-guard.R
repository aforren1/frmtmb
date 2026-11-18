# Review of lane slopetrig, check 6: the 5058 rule with its condition
# absent, beyond test-review-v29.R. True boundary correlations:
#   d2_novar   gaussian (1 + x | g), no slope variance (v29's design), n 200
#   d2_prop    gaussian (1 + x | g), slope effect = 0.6 * intercept effect
#              (true correlation 1), n 300
#   d3_prop    gaussian (1 + x1 + x2 | g), x2's effect = -0.5 * intercept's
#              (a Cholesky row at its limit), n 300
#   pois_prop  poisson (1 + x | g), slope effect = 0.5 * intercept's, n 600
#   binom_prop binomial(10) (1 + x | g), slope effect = 0.5 * intercept's
#   big_prop   gaussian (1 + x | g), 200 groups of 100 (20000 rows),
#              slope effect = 0.6 * intercept's: is a "gain" ever noise?
#   bigp_prop  poisson, 200 groups of 50 (10000 rows), as pois_prop
# Seeds 1..S, each at y * (1 + 2e-15 k), k = 0..K (counts: x moved).
# Per fit: correlation parameter past the clamp, gain_back, kind and lost
# reason, boundary message, "stopped short" text, and the shortfall
# against the best of D, T and a refit started with the correlation at 0.
# Usage: Rscript slopetrig-review-guard.R [S] [K]
a <- commandArgs(TRUE)
S <- if (length(a) > 0) as.integer(a[1]) else 5L
K <- if (length(a) > 1) as.integer(a[2]) else 2L
.libPaths(c("C:/Users/adf44/source/r/wt-slopetrig-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
mkd <- function(design, seed, k) {
  set.seed(seed)
  big <- design %in% c("big_prop", "bigp_prop")
  ng <- if (design == "big_prop") 200 else if (design == "bigp_prop") 200 else 20
  per <- if (design == "big_prop") 100 else if (design == "bigp_prop") 50 else
    if (design == "pois_prop") 30 else 15
  if (design == "d2_novar") per <- 10
  n <- ng * per
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(n)
  x2 <- rnorm(n)
  b0 <- rnorm(ng, 0, 0.7)
  eps <- 1 + 2e-15 * k
  out <- switch(design,
    d2_novar = data.frame(g = g, x = x,
                          y = (1 + 0.5 * x + b0[g] + rnorm(n)) * eps),
    d2_prop = , big_prop = data.frame(g = g, x = x,
                          y = (1 + 0.5 * x + b0[g] + 0.6 * b0[g] * x +
                                 rnorm(n)) * eps),
    d3_prop = data.frame(g = g, x1 = x, x2 = x2,
                         y = (1 + 0.5 * x + b0[g] + 0.3 * rnorm(ng)[g] * x -
                                0.5 * b0[g] * x2 + rnorm(n)) * eps),
    pois_prop = , bigp_prop = data.frame(g = g, x = x * eps,
                           y = rpois(n, exp(0.5 + 0.3 * x + 0.5 * b0[g] +
                                              0.25 * b0[g] * x))),
    binom_prop = data.frame(g = g, x = x * eps, n = 10L,
                            y = rbinom(n, 10, plogis(0.3 * x + b0[g] +
                                                       0.5 * b0[g] * x))))
  out
}
spec <- list(
  d2_novar = list(y ~ x + (1 + x | g), gaussian()),
  d2_prop = list(y ~ x + (1 + x | g), gaussian()),
  d3_prop = list(y ~ x1 + x2 + (1 + x1 + x2 | g), gaussian()),
  pois_prop = list(y ~ x + (1 + x | g), poisson()),
  binom_prop = list(y | trials(n) ~ x + (1 + x | g), binomial()),
  big_prop = list(y ~ x + (1 + x | g), gaussian()),
  bigp_prop = list(y ~ x + (1 + x | g), poisson()))
grab <- function(expr) {
  w <- character(0); m <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) {
      m <<- c(m, paste(class(x)[1], conditionMessage(x)))
      invokeRestart("muffleMessage")
    })
  list(v = v, w = w, m = m)
}
tc <- ns$se_cor_clamp_theta()
for (ds in names(spec)) {
  sp <- spec[[ds]]
  for (seed in seq_len(S)) for (k in 0:K) {
    d <- mkd(ds, seed, k)
    fd <- grab(frm(sp[[1]], data = d, family = sp[[2]]))
    if (inherits(fd$v, "error")) {
      cat(ds, seed, k, "ERROR", conditionMessage(fd$v), "\n"); next
    }
    f <- fd$v
    p <- f$opt$par
    thpos <- which(names(p) == "theta")
    bk <- f$frame$re_blocks[[1]]
    D <- bk$dim
    cpos <- thpos[-seq_len(D)]
    past <- abs(p[cpos]) > tc
    gains <- vapply(seq_along(cpos), function(i) {
      if (!past[i]) return(NA_real_)
      # the lane's own call: d = 2 one entry; wider, the Cholesky row
      if (D == 2) return(ns$se_cor_gain_back(f, cpos[i]))
      lt <- matrix(0L, D, D); lt[lower.tri(lt)] <- seq_along(cpos)
      r <- which(lt == i, arr.ind = TRUE)[1, 1]
      ns$se_cor_gain_back(f, cpos[lt[r, seq_len(r - 1L)]])
    }, 0)
    lost <- tryCatch(ns$sdr_of(f)$se_lost, error = function(e) NULL)
    kind <- f$cache$se_boundary_kind
    # a refit from the correlations at 0 (start = the fit's template)
    st <- par_template(f)
    st$theta[seq(D + 1, length(st$theta))] <- 0
    f0 <- grab(frm(sp[[1]], data = d, family = sp[[2]], start = unclass(st)))
    ft <- grab(frm(sp[[1]], data = d, family = sp[[2]],
                   control = frmtmb_control(autoscale = TRUE)))
    lls <- c(D = as.numeric(logLik(f)),
             S0 = if (inherits(f0$v, "error")) NA else as.numeric(logLik(f0$v)),
             T = if (inherits(ft$v, "error")) NA else as.numeric(logLik(ft$v)))
    cat(sprintf("%-10s s%d k%d past %s gain %s kind %s lost %s bmsg %d short %d nwarn %d shortfall %.3g\n",
                ds, seed, k, paste(past, collapse = ","),
                paste(signif(gains, 3), collapse = ","),
                paste(kind, collapse = ","),
                paste(paste0(names(lost), "=", lost), collapse = ","),
                sum(grepl("frmtmb_boundary_fit", fd$m)),
                sum(grepl("stopped short", c(fd$w, fd$m))), length(fd$w),
                max(lls, na.rm = TRUE) - lls[["D"]]))
  }
}
cat("DONE\n")
