# Punch round 1, minor 1: does engaging autoscale ever LOWER a
# skew_normal fit, and does the stationary-point escape take part?
# Default (engaged below sd 1e-3) against autoscale = FALSE, 12 seeds x
# spreads just under the threshold.
#   PREDFIX_ARM=lane Rscript dev/predfix-p1-skewdown.R > dev/predfix-log/p1-skewdown.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
mk <- function(seed, s, n = 250) {
  set.seed(seed)
  xs <- -abs(rnorm(n)) * 3
  yy <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = yy, xs = xs * s, xr = xs)
}
fos <- list(alpha0 = bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs),
            alpha1 = bf(y ~ xr, sigma ~ 1, alpha ~ xs),
            mu = bf(y ~ xs, sigma ~ 1, alpha ~ 1))
esc <- function(f) !is.null(f$opt[["stationary_escape"]])
cat(sprintf("%-7s %-4s %-8s %16s %16s %12s %5s %5s %6s %6s\n", "design", "seed", "s",
            "default", "FALSE", "def - FALSE", "cdef", "cF", "escD", "escF"))
rows <- list()
for (nm in names(fos)) for (seed in 1:12) for (s in c(5.4e-4, 2e-4, 1e-6)) {
  fo <- fos[[nm]]
  d <- mk(seed, s)
  a <- suppressWarnings(frm(fo, family = skew_normal(), data = d))
  b <- suppressWarnings(frm(fo, family = skew_normal(), data = d,
                            control = frmtmb_control(autoscale = FALSE)))
  dl <- as.numeric(logLik(a)) - as.numeric(logLik(b))
  rows[[length(rows) + 1L]] <- data.frame(design = nm, seed = seed, s = s, d = dl,
    ca = a$opt$convergence, cb = b$opt$convergence, ea = esc(a),
    eb = esc(b), engaged = !is.null(a$par_units))
  cat(sprintf("%-7s %-4d %-8g %16.9f %16.9f %12.3e %5d %5d %6s %6s\n", nm, seed, s,
              as.numeric(logLik(a)), as.numeric(logLik(b)), dl,
              a$opt$convergence, b$opt$convergence, esc(a), esc(b)))
}
r <- do.call(rbind, rows)
print(aggregate(cbind(lower = d < -1e-6, higher = d > 1e-6, engaged,
                      code_up = cb == 0 & ca != 0, escD = ea, escF = eb) ~
                  design, r, sum))
cat("\nengaged:", sum(r$engaged), "of", nrow(r),
    "; default lower than FALSE by > 1e-6:", sum(r$d < -1e-6),
    "; higher by > 1e-6:", sum(r$d > 1e-6),
    "; code 0 -> nonzero:", sum(r$cb == 0 & r$ca != 0),
    "; nonzero -> 0:", sum(r$cb != 0 & r$ca == 0),
    "; escape fired default/FALSE:", sum(r$ea), "/", sum(r$eb), "\n")
if (any(r$d < -1e-6)) print(r[r$d < -1e-6, ])
