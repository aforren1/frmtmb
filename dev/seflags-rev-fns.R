hyp_se <- function(fit, h) {
  r <- cap(hypothesis(fit, h, class = NULL))
  if (inherits(r$v, "rev_err")) return(sprintf("ERROR %s", r$v))
  sprintf("Est.Error %-12s warnings %d", format(r$v$hypothesis$Est.Error,
                                                digits = 6), length(r$w))
}

# jc_nonest() verdicts and propagated SEs for prediction gradients built
# at a chosen cosine with the stored flat direction (unit-diagonal
# metric), and their cosine with the EXACT null direction `nexact`
# (outer positions, natural units)
cos_sweep <- function(fit, nexact, label) {
  sdr <- suppressWarnings(ns$sdr_of(fit))
  fl <- sdr$se_null_flat
  on <- ns$outer_par_names(fit)
  cat(sprintf("  [%s] lost: %s\n", label, lost_of(fit)))
  if (is.null(fl)) {
    cat("  no flat subspace stored (se_null_flat NULL); null cols",
        NCOL(sdr$se_null), "\n")
    N <- sdr$se_null
    if (!NCOL(N)) return(invisible())
    v <- N[, NCOL(N)]
    sc <- sqrt(abs(diag(ns$fit_outer_hessian(fit)$H)))
    tol <- NA
  } else {
    N <- sdr$se_null
    v <- N[, fl$cols[1]]
    sc <- fl$scale
    tol <- fl$tol
  }
  r <- which(is.finite(sc) & sc > 0)
  vs <- v[r] * sc[r]
  vs <- vs / sqrt(sum(vs^2))
  ns_ <- nexact[r] * sc[r]
  ns_ <- ns_ / sqrt(sum(ns_^2))
  cat(sprintf("  flat tol %s | cosine(stored, exact) %.6f\n",
              format(tol, digits = 4), abs(sum(vs * ns_))))
  print(signif(stats::setNames(vs, on[r])[abs(vs) > 1e-3], 4))
  # w: orthogonal to the stored direction, within the parameters it loads
  ld <- abs(vs) > 0.05
  set.seed(11)
  w <- numeric(length(vs))
  w[ld] <- stats::rnorm(sum(ld))
  w <- w - sum(w * vs) * vs
  w <- w / sqrt(sum(w^2))
  pc <- ns$hyp_par_cov(fit)
  for (cc in c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8)) {
    gs <- cc * vs + sqrt(1 - cc^2) * w
    G <- numeric(length(on))
    G[r] <- gs * sc[r]
    gtrue <- abs(sum(gs * ns_)) / sqrt(sum(gs^2))
    bad <- ns$jc_nonest(list(null = N, flat = fl, units = fit$par_units),
                        matrix(G, 1), seq_along(on))
    vv <- suppressWarnings(ns$hyp_prop_var(pc, matrix(G[pc$outer_pos], 1)))
    cat(sprintf(paste0("  cos(stored) %.2f cos(exact) %.4f -> nonest %s",
                       " | hyp_prop_var se %s\n"),
                cc, gtrue, bad, format(sqrt(vv[1]), digits = 5)))
  }
}

