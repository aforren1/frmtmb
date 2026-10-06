# Reviewer, final check: on the RB4 fixture (dev/nanse-rev3-falseloss.R,
# seeds 11 and 17), compare se_tier3()'s global noise threshold with a
# per-direction one, computed here off the same Hessian and noise: for
# eigenvector v_k of the unit-diagonal free block, noise_k = ||Es v_k||.
# Which directions are "flat" under each rule, and which parameters
# would they take? (The lane's code is not changed; this reads its
# internals.)
#   Rscript dev/nanse-rev3-rb4-rule.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (s in c(11L, 17L)) {
  set.seed(s)
  n <- 480
  d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                  g2 = factor(rep(1:6, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$f] + rnorm(n)
  fit <- suppressWarnings(suppressMessages(
    frm(y ~ x + f + (1 + x | g2), data = d)))
  nm <- ns$outer_par_names(fit)
  h <- ns$fit_outer_hessian(fit)
  H <- h$H; E <- h$E
  u <- fit$par_units %||% rep(1, length(nm))
  Hu <- abs(H * outer(u, u)); Eu <- E * outer(u, u)
  empty <- apply(Hu, 1, max) <= pmax(10 * apply(Eu, 1, max), 1e-6)
  free <- which(!empty)
  Hf <- H[free, free]; Df <- sqrt(abs(diag(Hf)))
  e <- eigen(Hf / outer(Df, Df), symmetric = TRUE)
  Es <- E[free, free] / outer(Df, Df)
  big <- max(abs(e$values))
  glob <- max(1e-9 * big, 10 * sqrt(sum(Es^2)))
  perdir <- vapply(seq_along(e$values), function(k) {
    max(1e-9 * big, 10 * sqrt(sum((Es %*% e$vectors[, k])^2)))
  }, 0)
  cat(sprintf("seed %d: emptied %s | global flat_thr %.3g -> %d of %d directions removed\n",
              s, paste(nm[empty], collapse = ","), glob,
              sum(e$values <= glob), length(e$values)))
  rem <- which(e$values <= perdir)
  k0 <- which.min(e$values)
  cat(sprintf("   smallest eigenvalue %.3g, its own noise ||Es v|| %.3g (x10 = %.3g); free block min eigen > 1e-9 * largest: %s
",
              e$values[k0], perdir[k0] / 10, perdir[k0], min(e$values) > 1e-9 * big))
  cat(sprintf("   per-direction rule: %d removed (eigenvalues %s); dominant loadings %s\n",
              length(rem), paste(signif(e$values[rem], 3), collapse = " "),
              paste(unique(unlist(lapply(rem, function(k) {
                a <- abs(e$vectors[, k]); nm[free][a >= 0.5 * max(a)]
              }))), collapse = ",")))
}
