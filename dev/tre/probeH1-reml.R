# Probe H. Is REML defensible over a t latent?
#
# REML here means what it means everywhere in frmtmb: the mu fixed
# effects join the inner Laplace problem, so the objective is
#
#   L_R(sigma, s) = int prod_j m_j(beta) dbeta,
#   m_j(beta)     = int p(y_j | beta, b) p_t(b; s, nu) db,
#
# approximated by ONE Laplace step over (beta, b) jointly. The inner
# integral over b is the one probe A measured; the outer one over beta
# is new, and the question is whether stacking them degrades further.
#
# An exact reference exists for this design: `m_j` is the scalar
# quadrature A1 already certified, and beta is only two-dimensional, so
# the outer integral is a tensor adaptive Gauss-Hermite rule over the
# profile of the exact marginal likelihood.
#
# Run: Rscript dev/tre/probeH1-reml.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeH1.txt")
t0 <- Sys.time()

cat("Probe H: REML with t-distributed latents\n")
cat("========================================\n\n")

#' Exact REML log-likelihood: the beta integral taken by an adaptive
#' tensor Gauss-Hermite rule over the EXACT marginal likelihood.
exact_reml <- function(dat, sigma, s, nu, K = 25, Kin = 61) {
  ll <- function(beta) aghq_loglik(dat, beta, sigma, s, nu, K = Kin)
  op <- stats::optim(c(1, 0.5), function(b) -ll(b), method = "BFGS",
                     control = list(reltol = 1e-12))
  bhat <- op$par
  H <- numDeriv::hessian(function(b) -ll(b), bhat)
  R <- chol(H)                            # H = R'R
  gh <- statmod::gauss.quad(K, kind = "hermite")
  nodes <- as.matrix(expand.grid(gh$nodes, gh$nodes))
  lw <- log(gh$weights)
  lwt <- outer(lw, lw, "+")[cbind(rep(seq_len(K), K),
                                  rep(seq_len(K), each = K))]
  # z -> beta = bhat + sqrt(2) R^-1 z, jacobian 2^(d/2)/|R|
  B <- t(bhat + sqrt(2) * backsolve(R, t(nodes)))
  lv <- apply(B, 1, ll) + lwt + rowSums(nodes^2)
  mx <- max(lv)
  mx + log(sum(exp(lv - mx))) + log(2) - sum(log(diag(R)))
}

#' TMB REML: beta joins the random block, exactly as frmtmb does it.
fit_reml <- function(dat, nu, start = NULL, only_value = NULL) {
  ddat <- list(y = dat$y, X = dat$X, g = dat$g, nu_fix = nu)
  par <- list(log_sigma = 0, log_s = 0, beta = c(0, 0),
              b = rep(0, dat$G))
  fn <- function(p) {
    getAll(ddat, p)
    s <- exp(log_s)
    -sum(dnorm(y, X %*% beta + b[g], exp(log_sigma), log = TRUE)) -
      sum(dt(b / s, nu_fix, log = TRUE) - log_s)
  }
  obj <- MakeADFun(fn, par, random = c("beta", "b"), silent = TRUE)
  if (!is.null(only_value)) return(-obj$fn(only_value))
  op <- stats::nlminb(obj$par, obj$fn, obj$gr,
                      control = list(eval.max = 1000, iter.max = 500))
  list(par = c(sigma = exp(op$par[[1]]), s = exp(op$par[[2]])),
       loglik = -op$objective, obj = obj, opt = op)
}

fit_reml_gauss <- function(dat) {
  ddat <- list(y = dat$y, X = dat$X, g = dat$g)
  par <- list(log_sigma = 0, log_s = 0, beta = c(0, 0),
              b = rep(0, dat$G))
  fn <- function(p) {
    getAll(ddat, p)
    -sum(dnorm(y, X %*% beta + b[g], exp(log_sigma), log = TRUE)) -
      sum(dnorm(b, 0, exp(log_s), log = TRUE))
  }
  obj <- MakeADFun(fn, par, random = c("beta", "b"), silent = TRUE)
  op <- stats::nlminb(obj$par, obj$fn, obj$gr,
                      control = list(eval.max = 1000, iter.max = 500))
  list(par = c(sigma = exp(op$par[[1]]), s = exp(op$par[[2]])),
       loglik = -op$objective)
}

cat("A. REML objective at a fixed (sigma, s): TMB against the exact\n")
cat("   two-stage quadrature. G = 15 groups of n = 4.\n\n")
cat(sprintf("%-6s %14s %14s %12s %12s\n", "nu", "exact REML",
            "TMB REML", "error", "ML error"))
for (nu in c(2.5, 3, 5, 10, 30)) {
  d <- sim_tre(G = 15, n = 4, nu = nu, seed = 5)
  ex <- exact_reml(d, 1, 1, nu)
  tv <- fit_reml(d, nu, only_value = c(log_sigma = 0, log_s = 0))
  # the ML error on the same data, for scale
  mlex <- aghq_loglik(d, c(1, 0.5), 1, 1, nu, K = 101)
  mlla <- laplace_loglik(d, c(1, 0.5), 1, 1, nu)
  cat(sprintf("%-6s %14.6f %14.6f %12.5f %12.5f\n", nu, ex, tv,
              tv - ex, mlla - mlex))
}

cat("\nB. The whole REML fit: TMB against the exact objective\n")
cat("   maximized by grid refinement over (log sigma, log s).\n\n")
cat(sprintf("%-6s %-10s %12s %12s %12s\n", "nu", "fit", "sigma",
            "scale s", "logLik"))
for (nu in c(3, 5, 10)) {
  d <- sim_tre(G = 15, n = 4, nu = nu, seed = 5)
  ft <- fit_reml(d, nu)
  op <- stats::optim(log(ft$par), function(p)
    -exact_reml(d, exp(p[1]), exp(p[2]), nu, K = 21, Kin = 51),
    method = "Nelder-Mead",
    control = list(reltol = 1e-10, maxit = 400))
  cat(sprintf("%-6s %-10s %12.6f %12.6f %12.6f\n", nu, "exact",
              exp(op$par[1]), exp(op$par[2]), -op$value))
  cat(sprintf("%-6s %-10s %12.6f %12.6f %12.6f\n", nu, "TMB REML",
              ft$par[["sigma"]], ft$par[["s"]], ft$loglik))
  cat(sprintf("%17s difference: sigma %+.6f  scale %+.6f\n", "",
              ft$par[["sigma"]] - exp(op$par[1]),
              ft$par[["s"]] - exp(op$par[2])))
}

cat("\nC. The gaussian limit: REML with nu large against gaussian REML.\n\n")
cat(sprintf("%-12s %14s %14s %12s %12s\n", "nu", "sigma", "scale s",
            "d sigma", "d s"))
d <- sim_tre(G = 30, n = 6, nu = Inf, seed = 9)
fg <- fit_reml_gauss(d)
cat(sprintf("%-12s %14.9f %14.9f %12s %12s\n", "gaussian",
            fg$par[["sigma"]], fg$par[["s"]], "-", "-"))
for (nu in c(1e4, 1e6, 1e8)) {
  ft <- fit_reml(d, nu)
  cat(sprintf("%-12.0e %14.9f %14.9f %12.2e %12.2e\n", nu,
              ft$par[["sigma"]], ft$par[["s"]],
              abs(ft$par[["sigma"]] - fg$par[["sigma"]]),
              abs(ft$par[["s"]] - fg$par[["s"]])))
}

cat("\nD. Does REML move the scale in the right direction? Small G,\n")
cat("   where the ML variance component is most downward-biased.\n")
cat("   Mean over 60 replicates, true scale 1.\n\n")
cat(sprintf("%-6s %-6s %12s %12s\n", "nu", "G", "ML scale",
            "REML scale"))
for (nu in c(3, 5)) {
  for (G in c(8, 15, 40)) {
    ml <- reml <- numeric(60)
    for (r in 1:60) {
      dd <- sim_tre(G = G, n = 4, nu = nu, seed = 6000 + r)
      ml[r] <- fit_laplace(dd, nu = nu)$par[["s"]]
      reml[r] <- fit_reml(dd, nu)$par[["s"]]
    }
    cat(sprintf("%-6s %-6d %12.5f %12.5f\n", nu, G, mean(ml),
                mean(reml)))
  }
}

cat("\nelapsed: ", format(Sys.time() - t0), "\n")
sink()
cat("done\n")
