# Probe B. The inner Newton solve over a non-log-concave latent.
#
# A t log-density is CONCAVE only inside |b| < s*sqrt(nu); outside it is
# convex (its second derivative is +(nu+1)(b^2 - s^2 nu)/(s^2 nu + b^2)^2
# there, with a maximum of (nu+1)/(8 s^2 nu) at b^2 = 3 s^2 nu). The
# per-group conditional is that convexity plus the data's concavity
# -n/sigma^2, so a SECOND mode can only appear when
#
#     s^2 < sigma^2 (nu + 1) / (8 n nu),
#
# i.e. a latent scale small against the per-group residual precision.
# That is the low-information corner: singleton or near-singleton
# groups whose variance component is far below the residual. It is
# exactly where an outlying group's data pull one way and the t's spike
# at zero pulls the other. This probe (1) maps the region numerically,
# (2) checks whether TMB's inner optimizer finds the global mode there,
# and (3) measures what the Laplace error becomes when it does not.
#
# Run: Rscript dev/tre/probeB1-inner-newton.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeB1.txt")

cat("Probe B: multimodal per-group conditionals and the inner solve\n")
cat("==============================================================\n\n")

# --- B1: where is the conditional bimodal? ---------------------------
# Scan the per-group conditional h(b) on a dense grid and count local
# maxima directly, rather than trusting the analytic bound.
count_modes <- function(n, S, sigma, s, nu, lo = -20, hi = 20,
                        ngrid = 400001) {
  b <- seq(lo, hi, length.out = ngrid)
  st <- list(n = n, S = S, SS = S^2 / n)
  h <- h_mat(cbind(b), st, sigma, s, nu)[, 1]
  d <- diff(h)
  pk <- which(d[-length(d)] > 0 & d[-1] <= 0)
  list(k = length(pk), at = b[pk + 1], h = h[pk + 1])
}

cat("B1. Local maxima of the per-group conditional.\n")
cat("    sigma = 1, one group, residual mean rbar = S/n.\n")
cat("    'bound' is the analytic s below which two modes are possible.\n\n")
cat(sprintf("%-4s %-4s %-6s %-6s %8s %6s %-28s\n",
            "nu", "n", "s", "rbar", "bound", "modes", "mode locations"))
b1 <- NULL
for (nu in c(3, 5, 10)) {
  for (n in c(1, 2, 5)) {
    bound <- sqrt((nu + 1) / (8 * n * nu))
    for (s in c(0.15, 0.3, 0.6, 1.0)) {
      for (rbar in c(1, 3, 6)) {
        cm <- count_modes(n, rbar * n, 1, s, nu)
        cat(sprintf("%-4s %-4d %-6.2f %-6.1f %8.3f %6d %-28s\n",
                    nu, n, s, rbar, bound, cm$k,
                    paste(sprintf("%.3f", cm$at), collapse = ", ")))
        b1 <- rbind(b1, data.frame(nu = nu, n = n, s = s, rbar = rbar,
                                   bound = bound, modes = cm$k))
      }
    }
  }
}
cat("\nbimodal cases: ", sum(b1$modes > 1), " of ", nrow(b1), "\n", sep = "")
cat("all bimodal cases satisfy s < bound: ",
    all(b1$s[b1$modes > 1] < b1$bound[b1$modes > 1]), "\n", sep = "")

# --- B2: does TMB's inner Newton find the global mode? ---------------
cat("\n\nB2. TMB's inner solve on a dataset built inside that region.\n")
cat("    G groups of n = 2, sigma = 1, latent scale s, nu = 3, and\n")
cat("    ONE group forced to a far outlying residual mean.\n\n")

build_conflict <- function(G = 30, n = 2, s = 0.3, sigma = 1, nu = 3,
                           outlier = 6, seed = 7) {
  set.seed(seed)
  b <- s * rt(G, df = nu)
  b[1] <- 0                        # the outlier is in the DATA, not b
  d <- sim_tre(G = G, n = n, beta = c(0, 0), sigma = sigma, s = s,
               nu = nu, seed = seed, b_override = b)
  d$y[d$g == 1] <- d$y[d$g == 1] + outlier
  d
}

cat(sprintf("%-8s %-6s %10s %10s %10s %12s %12s\n",
            "outlier", "s", "TMB b1", "global", "other", "d(logLik)",
            "Laplace err"))
b2 <- NULL
for (s in c(0.15, 0.3, 0.6)) {
  for (outlier in c(2, 4, 6, 10)) {
    d <- build_conflict(s = s, outlier = outlier)
    fl <- tryCatch(fit_laplace(d, nu = 3), error = function(e) e)
    if (inherits(fl, "condition")) {
      cat(sprintf("%-8.1f %-6.2f  FIT ERROR: %s\n", outlier, s,
                  conditionMessage(fl)))
      next
    }
    p <- fl$par
    st <- suff_stats(d, c(p[["beta0"]], p[["beta1"]]))
    md <- h_mode(st, p[["sigma"]], p[["s"]], 3)
    # the inner solution TMB is using at its own optimum
    binner <- fl$b
    ex <- aghq_loglik(d, c(p[["beta0"]], p[["beta1"]]), p[["sigma"]],
                      p[["s"]], 3, K = 101, st = st)
    cat(sprintf("%-8.1f %-6.2f %10.4f %10.4f %10.4f %12.2e %12.4f\n",
                outlier, s, binner[1], md$mode[1], md$other[1],
                max(abs(binner - md$mode)), fl$loglik - ex))
    b2 <- rbind(b2, data.frame(
      s = s, outlier = outlier, tmb_b1 = binner[1],
      global_b1 = md$mode[1], other_b1 = md$other[1],
      maxdiff = max(abs(binner - md$mode)),
      n_bimodal = sum(md$bimodal), laplace_err = fl$loglik - ex))
  }
}

# --- B3: force the wrong mode and price it ---------------------------
cat("\n\nB3. The price of landing on the wrong mode.\n")
cat("    A single group, both modes taken in turn.\n\n")
cat(sprintf("%-4s %-4s %-6s %-6s %10s %10s %10s %10s\n",
            "nu", "n", "s", "rbar", "mode A", "mode B",
            "Lap at A", "Lap at B"))
for (nu in c(3, 5)) {
  for (s in c(0.15, 0.25)) {
    for (rbar in c(2, 4, 8)) {
      n <- 1
      cm <- count_modes(n, rbar * n, 1, s, nu)
      if (cm$k < 2) next
      st <- list(n = n, S = rbar * n, SS = (rbar * n)^2 / n)
      lap_at <- function(bb) {
        hh <- h_mat(cbind(bb), st, 1, s, nu)[1, 1]
        cv <- -( -n / 1 + ldt_hess(bb, s, nu))
        hh + 0.5 * log(2 * pi) - 0.5 * log(cv)
      }
      cat(sprintf("%-4s %-4d %-6.2f %-6.1f %10.4f %10.4f %10.4f %10.4f\n",
                  nu, n, s, rbar, cm$at[1], cm$at[2],
                  lap_at(cm$at[1]), lap_at(cm$at[2])))
      cat(sprintf("%36s exact %10.4f\n", "",
                  log(stats::integrate(function(b)
                    exp(h_mat(cbind(b), st, 1, s, nu)[, 1]),
                    -Inf, Inf, rel.tol = 1e-11)$value)))
    }
  }
}

# --- B4: non-positive-definite inner Hessian ------------------------
cat("\n\nB4. Does the inner Hessian ever leave the PD cone during the\n")
cat("    outer optimization?  Sweep log_s over a wide range at the\n")
cat("    conflict data and record TMB's inner behaviour.\n\n")
d <- build_conflict(s = 0.3, outlier = 8)
fl <- fit_laplace(d, nu = 3, silent = TRUE)
p <- fl$opt$par
cat(sprintf("%-10s %14s %14s %10s\n", "log_s", "TMB -logLik",
            "hand Laplace", "min curv"))
for (ls in seq(-4, 2, by = 0.5)) {
  pp <- p
  pp[["log_s"]] <- ls
  v <- tryCatch(fl$obj$fn(pp), error = function(e) NA_real_)
  st <- suff_stats(d, pp[1:2])
  md <- h_mode(st, exp(pp[["log_sigma"]]), exp(ls), 3)
  hand <- laplace_loglik(d, pp[1:2], exp(pp[["log_sigma"]]), exp(ls), 3,
                         st = st)
  cat(sprintf("%-10.2f %14.6f %14.6f %10.4f\n", ls, -v, hand,
              min(md$curv)))
}

saveRDS(list(b1 = b1, b2 = b2), "dev/tre/probeB1.rds")
sink()
cat("done\n")
