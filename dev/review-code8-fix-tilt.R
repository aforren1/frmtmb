# Fix-round check (B1): try to break the ritz_flat exception. The ridge
# c0 + exp(a)^k made gently tilted by a prior on k (MAP objective), so
# the direction the Ritz step calls flat is not flat: the stop the
# optimizer leaves on the ridge can be short along it. Reference: the
# best of nlminb and BFGS polishes started at the stop and at points
# along the ridge with k = 0.2, 0.5, 1.
#   RV_ARM=before|after Rscript dev/review-code8-fix-tilt.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
rows <- list()
for (seed in c(77, 1, 2, 3)) {
  set.seed(seed)
  dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(8, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  y0 <- dn$yn
  for (psd in c(30, 10, 3, 1)) for (j in c(0, 5, -9)) {
    dn$yn <- y0 * (1 + 2e-15 * j)
    pr <- set_prior(sprintf("normal(0.5, %g)", psd), nlpar = "k")
    r <- rfit(form, data = dn, prior = pr)
    if (r$err) next
    f <- r$f; obj <- f$obj; x <- f$opt$par
    nm <- ns$outer_par_names(f)
    ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
    f0 <- obj$fn(x)
    best <- f0
    starts <- list(x)
    for (kk in c(0.2, 0.5, 1)) {
      z <- x; z[ik] <- kk; z[ia] <- x[ia] * x[ik] / kk
      starts[[length(starts) + 1L]] <- z
    }
    for (z in starts) {
      o1 <- tryCatch(nlminb(z, obj$fn, obj$gr,
                            control = list(iter.max = 3000, eval.max = 6000)),
                     error = function(e) list(objective = Inf, par = z))
      o2 <- tryCatch(optim(o1$par, obj$fn, obj$gr, method = "BFGS",
                           control = list(maxit = 5000, reltol = 1e-14)),
                     error = function(e) list(value = Inf))
      best <- min(best, o1$objective, o2$value, na.rm = TRUE)
    }
    s <- se_of(f)
    rows[[length(rows) + 1L]] <- data.frame(seed = seed, psd = psd, j = j,
      code = f$opt$convergence, k = signif(x[ik], 3),
      short = signif(f0 - best, 3),
      grad_w = any(grepl("Large maximum", r$w)),
      code_w = any(r$wc == "frmtmb_conv_code"),
      at_max_m = any(r$mc == "frmtmb_converged_at_max"),
      lost = paste(sort(names(s$lost)), collapse = ","),
      nwarn = length(r$w))
  }
}
X <- do.call(rbind, rows)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/fix-tilt-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
bad <- X$short > 1e-3 & !X$grad_w & !X$code_w
cat("short > 1e-3:", sum(X$short > 1e-3), " with no gradient/code warning:",
    sum(bad), " shortfalls:", paste(X$short[bad], collapse = " "), "\n")
cat("DONE\n")
