# Lane nlstop, item 4705: the rows of dev/review-code8-fix-ridgeprior.R
# and -ridgeprior2.R on one table. The ridge c0 + exp(a)^k with
# normal(0.5, psd) on k is tilted along the ridge by the prior only. A
# custom optimizer runs nlminb to the MAP, then moves exactly along the
# ridge (a, k, theta_1) -> (s a, k / s, theta_1 + log s) to k = kt and
# reports code 8 (1) or 0. A stop more than grad_tol (1e-3) short must
# warn.
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-ridgeprior.R
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
grid <- list(`10` = c(5, 20, 50, 200), `100` = c(5, 20, 50, 200),
             `300` = c(50, 80, 110, 150),
             `1000` = c(5, 20, 50, 80, 110, 150, 200))
cat("ROWHEAD psd code kt short grad_w at_max_m says_max nwarn classes\n")
for (psd in names(grid)) {
  pr <- set_prior(sprintf("normal(0.5, %s)", psd), nlpar = "k")
  f <- suppressWarnings(suppressMessages(frm(form, data = dn, prior = pr)))
  nm <- ns$outer_par_names(f)
  ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
  it <- match("theta_1", nm)
  for (code in c(1L, 0L)) for (kt in grid[[psd]]) {
    opt <- function(par, fn, gr, lower, upper, control) {
      o <- stats::nlminb(par, fn, gr, control = list(iter.max = 1000,
                                                     eval.max = 2000))
      p <- o$par
      s <- p[ik] / kt
      p[ia] <- p[ia] * s; p[ik] <- kt; p[it] <- p[it] + log(abs(s))
      list(par = p, objective = fn(p), convergence = code,
           message = if (code == 1L) "false convergence (8)" else
             "relative convergence (4)")
    }
    r <- tt(frm(form, data = dn, prior = pr,
                control = frmtmb_control(optimizer = opt, restarts = 0)))
    if (iserr(r)) { cat("ERR", r$msg, "\n"); next }
    w <- attr(r, "w")
    cat(sprintf("ROW %s %d %g %.4g %s %s %s %d %s\n", psd, code, kt,
                r$opt$objective - f$opt$objective,
                any(grepl("Large maximum", w)),
                any(attr(r, "mc") == "frmtmb_converged_at_max"),
                any(grepl("so the estimates are at the maximum",
                          c(w, attr(r, "m")))),
                length(w), wsum(r)))
  }
}
cat("DONE\n")
