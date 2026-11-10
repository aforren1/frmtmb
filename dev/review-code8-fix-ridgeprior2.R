# Third fix-round check: the Ritz-flat exception under a weak prior. The
# ridge c0 + exp(a)^k with normal(0.5, psd) on k: the objective is tilted
# along the ridge by the prior only. An optimizer runs nlminb to the MAP
# and then moves ALONG the ridge (a, k, theta_1) -> (s a, k / s,
# theta_1 + log s), exactly, to k = k_target; reported code 8 or 0.
#   RV_ARM=before|after Rscript dev/review-code8-fix-ridgeprior.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
rows <- list()
for (psd in c(300, 1000)) {
  pr <- set_prior(sprintf("normal(0.5, %g)", psd), nlpar = "k")
  ref <- rfit(form, data = dn, prior = pr)
  f <- ref$f
  nm <- ns$outer_par_names(f)
  ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
  it <- match("theta_1", nm)
  for (code in c(1L, 0L)) for (kt in c(50, 80, 110, 150)) {
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
    r <- rfit(form, data = dn, prior = pr,
              control = frmtmb_control(optimizer = opt, restarts = 0))
    if (r$err) { cat("err", conditionMessage(r$f), "\n"); next }
    g <- r$f
    txt <- c(r$w, r$m)
    rows[[length(rows) + 1L]] <- data.frame(psd = psd, code = code, kt = kt,
      short = signif(g$opt$objective - f$opt$objective, 3),
      says_max = any(grepl("so the estimates are at the maximum", txt)),
      grad_w = any(grepl("Large maximum", r$w)),
      at_max_m = any(r$mc == "frmtmb_converged_at_max"),
      nwarn = length(r$w),
      lost = paste(sort(names(se_of(g)$lost)), collapse = ","))
  }
}
X <- do.call(rbind, rows)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/fix-ridgeprior2-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DONE\n")
