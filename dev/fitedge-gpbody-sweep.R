# Item 1 sweep: gp() inside a nonlinear body. 4 bodies x 4 curves x n 60, 120 x
# seeds 1 to 3 = 96 rows. Args: arm (lane | base) and part:
#   fit  the default fit of the body, one line per row
#   lin  the default fit of the linear spelling (bodies A and B only)
#   ref  the reference: the best of the linear spelling (where one exists) and
#        12 theta starts on the body (a theta start skips the warm start, so
#        the reference does not depend on the build)
# Seeds: set.seed(18000 + 100 * seed + n + curve index).
a <- commandArgs(TRUE)
arm <- a[1]; part <- a[2]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
curves <- c("sin", "wig", "exp", "step")
mk <- function(seed, n, ci, body) {
  set.seed(18000 + 100 * seed + n + ci)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
  f <- switch(curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  y <- switch(body,
    A = f + 0.3 * z + 0.1 * rnorm(n),
    B = rpois(n, exp(1 + 0.5 * f / max(abs(f)) + 0.3 * z)),
    C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n),
    D = 1 + f * w + 0.1 * rnorm(n))
  data.frame(x, z, w, y)
}
spell <- list(
  A = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = y ~ z + gp(x), fam = gaussian()),
  B = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = y ~ z + gp(x), fam = poisson()),
  C = list(nl = bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = NULL, fam = gaussian()),
  D = list(nl = bf(y ~ b + a * w, a ~ gp(x), b ~ 1, nl = TRUE),
           lin = NULL, fam = gaussian()))
out <- file.path(FE_OUT, sprintf("gpbody-%s-%s.tsv", part, arm))
cat("body\tcurve\tn\tseed\tll\tconv\tnwarn\twclass\tsec\n", file = out)
for (body in names(spell)) for (ci in seq_along(curves)) for (n in c(60L, 120L))
for (seed in 1:3) {
  d <- mk(seed, n, ci, body)
  sp <- spell[[body]]
  t0 <- proc.time()[[3]]
  if (part %in% c("fit", "lin")) {
    if (part == "lin" && is.null(sp$lin)) next
    r <- fe_fit(if (part == "lin") sp$lin else sp$nl, data = d, family = sp$fam)
    ll <- if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit))
    cv <- if (inherits(r$fit, "error")) NA else r$fit$opt$convergence
    wc <- if (inherits(r$fit, "error")) "ERROR" else
      paste(unique(sub(":.*", "", r$warn)), collapse = ",")
    nw <- length(r$warn)
  } else {
    best <- -Inf
    if (!is.null(sp$lin)) {
      r <- fe_fit(sp$lin, data = d, family = sp$fam)
      if (!inherits(r$fit, "error")) best <- as.numeric(logLik(r$fit))
    }
    for (rho in c(0.15, 0.5, 1.5, 5)) for (sd0 in c(0.3, 1, 3)) {
      th <- c(log(sd0), log(rho))
      r <- fe_fit(sp$nl, data = d, family = sp$fam, start = list(theta = th))
      if (!inherits(r$fit, "error")) {
        v <- as.numeric(logLik(r$fit))
        if (is.finite(v) && v > best) best <- v
      }
    }
    ll <- best; cv <- NA; nw <- NA; wc <- ""
  }
  cat(sprintf("%s\t%s\t%d\t%d\t%.6f\t%s\t%s\t%s\t%.2f\n", body, curves[ci], n,
              seed, ll, cv, nw, wc, proc.time()[[3]] - t0), file = out,
      append = TRUE)
}
cat("DONE", out, "\n")
