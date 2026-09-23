# Two questions about the skew_normal() stall, measured rather than
# assumed:
#
# 1. The rate, on the design exactly as dev/drmtmb-findings.md states
#    it ("stated" below) and on the stream the first run used, which
#    drew one extra rnorm(n) before the covariate ("with dead draw").
# 2. Whether either candidate fix removes the stall: start alpha from
#    the RESIDUAL skewness (frmtmb's own rule, applied to lm residuals
#    instead of the raw response), or refit from +2 and -2 and keep the
#    better optimum. Neither is implemented in R/; both are applied here
#    through frm(start = ).
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
skew <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
start_rule <- function(m3) 2 * sign(m3) + 0.5 * m3

make_data <- function(seed, dead_draw) {
  set.seed(seed)
  n <- 200
  if (dead_draw) x <- rnorm(n)      # the first run's stream
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = y, xs = xs)
}
# Only alpha moves. frmtmb's own sigma start is log(sd(y)), and passing
# a different one changes the outcome: at seed 34 of the dead-draw
# stream, the residual-skewness start reaches the optimum with sigma
# started at 0 and stalls with frmtmb's own sigma start. An in-family
# fix would change alpha alone, so that is what this measures.
fit_at <- function(dd, a0, sigma0 = log(stats::sd(dd$y))) {
  frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = dd,
      start = list(betad = c(sigma0, a0)))
}
one_seed <- function(seed, dead_draw) {
  dd <- make_data(seed, dead_draw)
  ref <- drmTMB::drmTMB(dbf(y ~ xs, sigma ~ 1, nu ~ 1),
                        family = drmTMB::skew_normal(), data = dd)
  ll_ref <- as.numeric(logLik(ref))
  def <- frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(),
             data = dd)
  a_res <- start_rule(skew(residuals(stats::lm(y ~ xs, data = dd))))
  res <- fit_at(dd, a_res)
  two <- list(fit_at(dd, 2), fit_at(dd, -2))
  ll_two <- max(vapply(two, function(f) as.numeric(logLik(f)), 0))
  c(seed = seed, ll_ref = ll_ref,
    gap_default = as.numeric(logLik(def)) - ll_ref,
    gap_resid_start = as.numeric(logLik(res)) - ll_ref,
    gap_two_start = ll_two - ll_ref,
    alpha_resid_start = a_res,
    alpha_resid_fit = res$opt$par[[4]])
}
report <- function(dead_draw) {
  lab <- if (dead_draw) "stream WITH the dead rnorm(n) draw" else
    "stream AS STATED in dev/drmtmb-findings.md"
  cat("\n==========", lab, "\n")
  r <- as.data.frame(t(vapply(1:40, one_seed, numeric(7),
                              dead_draw = dead_draw)))
  print(round(r, 4), row.names = FALSE)
  for (nm in c("gap_default", "gap_resid_start", "gap_two_start")) {
    bad <- r[[nm]] < -1e-6
    cat(sprintf("%-16s stalls: %2d of 40", nm, sum(bad)))
    if (any(bad)) {
      cat(sprintf("   gaps %.6f to %.6f   seeds %s",
                  -max(r[[nm]][bad]), -min(r[[nm]][bad]),
                  paste(r$seed[bad], collapse = " ")))
    }
    cat("\n")
  }
  invisible(r)
}
r_stated <- report(FALSE)
r_dead <- report(TRUE)

# The reviewer's seed 34 of the dead-draw stream: a correctly signed
# residual start that still falls into alpha = 0. Scan starts there.
s34 <- r_dead[r_dead$seed == 34, ]
cat("\nseed 34, dead-draw stream: residual start",
    format(s34$alpha_resid_start, digits = 6), "-> alpha",
    format(s34$alpha_resid_fit, digits = 6), " gap",
    format(s34$gap_resid_start, digits = 6), "\n")
dd34 <- make_data(34, TRUE)
cat("  alpha start / logLik with frmtmb's sigma start / with sigma 0\n")
for (a0 in c(0.5, 1, 2, 2.5524, 3, 4, 5, 8, 10, 16)) {
  f <- fit_at(dd34, a0)
  g <- fit_at(dd34, a0, sigma0 = 0)
  cat(sprintf("  %7.4f  %12.4f (alpha %9.5f)  %12.4f\n", a0,
              as.numeric(logLik(f)), f$opt$par[[4]],
              as.numeric(logLik(g))))
}
