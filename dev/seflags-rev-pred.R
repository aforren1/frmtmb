# Reviewer of lane seflags, claim 1: predictions along a lost direction
# must stay NaN, and a prediction with a moderate projection onto an
# EXACTLY flat direction is not estimable whatever the threshold says.
#   Rscript dev/seflags-rev-pred.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")

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

## A. Spread ridge, k = 10, exact Hessian (nanse B2)
spread <- function(k, m = 5, seed = 1, re = FALSE) {
  set.seed(seed)
  f <- factor(rep(seq_len(k), each = m))
  d <- data.frame(f = f, h = factor(rep(1:5, length.out = k * m)),
                  y = rnorm(k)[f] + rnorm(k * m, 0, 0.5))
  if (re) d$y <- d$y + rnorm(5, 0, 0.5)[d$h]
  d
}
for (re in c(FALSE, TRUE)) {
  dd <- spread(10, re = re)
  fo <- if (re) {
    bf(y ~ a + b, a ~ 0 + f, b ~ 1 + (1 | h), nl = TRUE)
  } else bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE)
  r <- cap(frm(fo, data = dd))
  fit <- r$v
  cat("\n== A. spread k = 10", if (re) "with b ~ 1 + (1 | h)" else "exact",
      "| code", fit$opt$convergence, "| fit warnings", length(r$w), "\n")
  nd <- dd[c(1, 6), , drop = FALSE]
  for (dp in list(NULL, "a", "b")) {
    p <- cap(fitted(fit, newdata = nd, dpar = dp, re_formula = NA))
    cat(sprintf("  fitted(dpar = %s): Est.Error %s warnings %d\n",
                if (is.null(dp)) "mu" else dp,
                paste(signif(p$v[, 2], 4), collapse = " "), length(p$w)))
  }
  sig <- exp(fit$opt$par[length(fit$opt$par)])
  for (h in c("b_a_f1 + b_b_Intercept = 0", "b_a_f1 = 0",
              "b_a_f1 - b_a_f2 = 0", "b_a_f1 + 0.5 * b_b_Intercept = 0",
              "b_a_f1 + b_a_f2 + b_b_Intercept = 0",
              "b_a_f1 + b_a_f2 + b_a_f3 + 2 * b_b_Intercept = 0",
              "b_a_f1 + b_a_f2 + 1.5 * b_b_Intercept = 0")) {
    cat(sprintf("  %-46s %s\n", h, hyp_se(fit, h)))
  }
  on <- ns$outer_par_names(fit)
  nex <- numeric(length(on))
  nex[grepl("^a_f", on)] <- 1
  nex[on == "b_(Intercept)"] <- -1
  cos_sweep(fit, nex, "spread")
}

## B. y ~ x + s(x) + (1 | g), the lane's alias construction, means 1e-5
## and 1e-3 (intercept kept / lost)
for (m in c(1e-5, 1e-3, 0.1)) {
  set.seed(1)
  n <- 200
  x <- stats::rnorm(n)
  x <- x - mean(x) + m
  d <- data.frame(x = x, g = factor(rep(1:10, n / 10)))
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + stats::rnorm(10, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.3)
  r <- cap(frm(y ~ x + s(x) + (1 | g), data = d))
  fit <- r$v
  cat("\n== B. alias, mean", m, "| code", fit$opt$convergence,
      "| fit warnings", length(r$w), "\n")
  X <- fit$frame$linpreds[[1]]$X
  cn <- colnames(X)
  ifx <- grep("^s[(]x[)][.]?fx", cn)
  cf <- coef(stats::lm(as.numeric(X[, ifx[1]]) ~ d$x))
  on <- ns$outer_par_names(fit)
  nex <- numeric(length(on))
  nex[on == "(Intercept)"] <- cf[1]
  nex[on == "x"] <- cf[2]
  nex[grepl("fx1$", on)] <- -1
  cat("  exact alias: fx1 =", format(cf, digits = 6), "\n")
  print(variables(fit)[1:6])
  cos_sweep(fit, nex, "alias")
  # user-facing: the two named coefficients mixed so that the cosine with
  # the alias is moderate
  vn <- variables(fit)
  bx <- vn[grepl("^b_x$", vn)]
  bs <- vn[grepl("^bs_sx_1$", vn)]
  if (length(bx) && length(bs)) {
    for (h in c(sprintf("%s = 0", bx),
                sprintf("%s + %g * %s = 0", bx, 1 / cf[2], bs),
                sprintf("%s + %g * %s = 0", bx, 0.5 / cf[2], bs),
                sprintf("%s + %g * %s = 0", bx, 0.8 / cf[2], bs),
                sprintf("%s + %g * %s = 0", bx, 0.9 / cf[2], bs))) {
      cat(sprintf("  %-46s %s\n", h, hyp_se(fit, h)))
    }
  }
}
