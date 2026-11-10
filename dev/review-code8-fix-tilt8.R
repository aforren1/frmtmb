# Fix-round check (B1), forced: the gently tilted ridge (prior normal(0.5,
# psd) on k), stopped early by optim()'s BFGS at maxit = k and reported
# as "false convergence (8)" (the device of conv_fake8()). A stop on the
# ridge short along its tilt is where the ritz_flat exception could fire.
#   RV_ARM=after Rscript dev/review-code8-fix-tilt8.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
mk8 <- function(k) {
  function(par, fn, gr, lower, upper, control) {
    o <- stats::optim(par, fn, gr, method = "BFGS", control = list(maxit = k))
    list(par = o$par, objective = o$value, convergence = 1L,
         message = "false convergence (8)")
  }
}
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
rows <- list()
for (seed in c(77, 1, 2)) {
  set.seed(seed)
  dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(8, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  for (psd in c(100, 30, 10, Inf)) {
    pr <- if (is.finite(psd)) {
      set_prior(sprintf("normal(0.5, %g)", psd), nlpar = "k")
    }
    ref <- rfit(form, data = dn, prior = pr)
    fref <- ref$f$opt$objective
    for (k in c(15, 20, 30, 40, 60, 80, 120, 200)) {
      r <- rfit(form, data = dn, prior = pr,
                control = frmtmb_control(optimizer = mk8(k), restarts = 0))
      if (r$err) next
      f <- r$f
      txt <- c(r$w, r$m)
      rows[[length(rows) + 1L]] <- data.frame(seed = seed, psd = psd, k = k,
        short = signif(f$opt$objective - fref, 3),
        says_max = any(grepl("so the estimates are at the maximum", txt)),
        grad_w = any(grepl("Large maximum", r$w)),
        at_max_m = any(r$mc == "frmtmb_converged_at_max"),
        nwarn = length(r$w),
        lost = paste(sort(names(se_of(f)$lost)), collapse = ","))
    }
  }
}
X <- do.call(rbind, rows)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/review-code8-out/fix-tilt8-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("short > 2e-5 with the message:", sum(X$at_max_m & X$short > 2e-5),
    "of", sum(X$short > 2e-5), "\n")
cat("short > 1e-3 called at max:", sum(X$says_max & X$short > 1e-3),
    "; short > 1e-3 without gradient warning:",
    sum(!X$grad_w & X$short > 1e-3), "\n")
cat("DONE\n")
