# Punch round 1, M1: gaussian y ~ x + (1 + x | g) with x rescaled. ML is
# invariant to the scale of x, so every scale should reach the scale-1
# optimum. Where does each arm stop, and what is theta there?
#   PREDFIX_ARM=lane Rscript dev/predfix-p1-slope.R > dev/predfix-log/p1-slope-<tag>.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
mk <- function(seed, n_g = 20, n_per = 15) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n_g * n_per)
  u0 <- rnorm(n_g, 0, 0.7)
  u1 <- rnorm(n_g, 0, 0.4)
  data.frame(g = g, x = x, y = 1 + 0.5 * x + u0[g] + u1[g] * x +
               rnorm(n_g * n_per))
}
fo <- bf(y ~ x + (1 + x | g))
ll <- function(f) if (inherits(f, "try-error")) NA else as.numeric(logLik(f))
cat(sprintf("%-4s %-8s %-7s %14s %12s %5s %8s %s\n", "seed", "sd(x)", "arm",
            "logLik", "short", "code", "engaged", "theta"))
for (seed in 1:3) {
  d0 <- mk(seed)
  ref <- ll(frm(fo, family = gaussian(), data = d0))
  for (s in c(1e-1, 3e-2, 1e-2, 1e-3, 1e-4, 1e-6)) {
    d <- d0
    d$x <- d0$x * s
    for (a in list(NULL, FALSE, TRUE)) {
      f <- try(suppressWarnings(frm(fo, family = gaussian(), data = d,
                                    control = frmtmb_control(autoscale = a))),
               silent = TRUE)
      th <- if (inherits(f, "try-error")) NA else f$estimates$theta
      cat(sprintf("%-4d %-8.3g %-7s %14.6f %12.6f %5s %8s %s\n", seed,
                  sd(d$x), format(if (is.null(a)) "NULL" else a), ll(f),
                  ref - ll(f),
                  if (inherits(f, "try-error")) "err" else f$opt$convergence,
                  if (inherits(f, "try-error")) NA else !is.null(f$par_units),
                  paste(signif(th, 4), collapse = " ")))
    }
  }
}
