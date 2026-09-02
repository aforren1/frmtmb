# Probe G. Can `quadrature = TRUE` marginalize a t latent EXACTLY?
#
# frmtmb's quadrature path is TMB's experimental `integrate =`
# Gauss-Kronrod marginalization over each scalar random effect. The rule
# is distribution-agnostic in principle, since it integrates whatever
# the tape produces. But the rescaling is calibrated ONCE at the Laplace
# mode and curvature and then frozen, which is the part a fat tail can
# defeat: a rule sized for a gaussian shoulder truncates a polynomial
# tail.
#
# This is worth settling before implementation, because if it holds it
# turns the whole Laplace-accuracy caveat into a one-argument check the
# user can run.
#
# Run: Rscript dev/tre/probeG1-quadrature.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeG1.txt")

gk <- list(b = structure(list(dim = 1, adaptive = FALSE, debug = FALSE,
                              method = "marginal_gk"), class = "GK"))

#' The same model as fit_laplace(), built with integrate = GK.
fit_gk <- function(dat, nu, anchor) {
  ddat <- list(y = dat$y, X = dat$X, g = dat$g, nu_fix = nu)
  fn <- function(p) {
    getAll(ddat, p)
    s <- exp(log_s)
    -sum(dnorm(y, X %*% beta + b[g], exp(log_sigma), log = TRUE)) -
      sum(dt(b / s, nu_fix, log = TRUE) - log_s)
  }
  obj <- MakeADFun(fn, anchor, random = "b", integrate = gk,
                   silent = TRUE)
  f0 <- try(obj$fn(obj$par), silent = TRUE)
  if (inherits(f0, "try-error") || !is.finite(f0)) return(NULL)
  op <- try(stats::nlminb(obj$par, obj$fn, obj$gr,
                          control = list(eval.max = 1000,
                                         iter.max = 500)),
            silent = TRUE)
  if (inherits(op, "try-error")) return(NULL)
  list(par = c(beta0 = op$par[[1]], beta1 = op$par[[2]],
               sigma = exp(op$par[[3]]), s = exp(op$par[[4]])),
       loglik = -op$objective, obj = obj, opt = op)
}

cat("Probe G: Gauss-Kronrod quadrature over a t latent\n")
cat("=================================================\n\n")

cat("A. logLik at a FIXED parameter vector (the truth), against the\n")
cat("   AGHQ reference and the Laplace value. 40 groups.\n\n")
cat(sprintf("%-6s %-4s %14s %14s %14s %12s %12s\n", "nu", "n",
            "exact (AGHQ)", "Laplace", "GK", "GK err", "Lap err"))
for (nu in c(2.5, 3, 5, 10)) {
  for (n in c(3, 8)) {
    d <- sim_tre(G = 40, n = n, nu = nu, seed = 11)
    st <- suff_stats(d, c(1, 0.5))
    ex <- aghq_loglik(d, c(1, 0.5), 1, 1, nu, K = 101, st = st)
    la <- laplace_loglik(d, c(1, 0.5), 1, 1, nu, st = st)
    # calibrate the GK tape at the conditional modes for THESE values
    fl <- fit_laplace(d, nu = nu)
    p <- c(1, 0.5, 0, 0)
    fl$obj$fn(p)                       # solves the inner problem there
    anchor <- fl$obj$env$parList(p, fl$obj$env$last.par)
    obj <- tryCatch({
      ddat <- list(y = d$y, X = d$X, g = d$g, nu_fix = nu)
      fn <- function(pp) {
        getAll(ddat, pp)
        s <- exp(log_s)
        -sum(dnorm(y, X %*% beta + b[g], exp(log_sigma), log = TRUE)) -
          sum(dt(b / s, nu_fix, log = TRUE) - log_s)
      }
      anchor$log_nu <- NULL
      MakeADFun(fn, anchor[c("beta", "log_sigma", "log_s", "b")],
                random = "b", integrate = gk, silent = TRUE)
    }, error = function(e) e)
    gkv <- if (inherits(obj, "condition")) NA_real_ else
      -obj$fn(c(1, 0.5, 0, 0))
    cat(sprintf("%-6s %-4d %14.6f %14.6f %14.6f %12.5f %12.5f\n",
                nu, n, ex, la, gkv, gkv - ex, la - ex))
  }
}

cat("\n\nB. Full fits: GK vs Laplace vs exact AGHQ ML, same data.\n\n")
cat(sprintf("%-6s %-4s %-10s %10s %10s %10s %12s\n", "nu", "n", "fit",
            "beta0", "sigma", "scale s", "logLik"))
for (nu in c(2.5, 3, 5, 10)) {
  for (n in c(3, 8)) {
    d <- sim_tre(G = 40, n = n, nu = nu, seed = 11)
    fl <- fit_laplace(d, nu = nu)
    fe <- fit_exact(d, nu = nu, K = 101,
                    start = c(fl$par[1], fl$par[2], log(fl$par[3]),
                              log(fl$par[4])))
    fl$obj$fn(fl$opt$par)
    anchor <- fl$obj$env$parList(fl$opt$par, fl$obj$env$last.par)
    anchor$log_nu <- NULL
    fg <- tryCatch(fit_gk(d, nu, anchor[c("beta", "log_sigma", "log_s",
                                          "b")]),
                   error = function(e) NULL)
    pr <- function(lab, p, ll) {
      cat(sprintf("%-6s %-4d %-10s %10.5f %10.5f %10.5f %12.5f\n",
                  nu, n, lab, p[[1]], p[[3]], p[[4]], ll))
    }
    pr("exact", fe$par, fe$loglik)
    pr("Laplace", fl$par, fl$loglik)
    if (is.null(fg)) {
      cat(sprintf("%-6s %-4d %-10s  GK FIT FAILED\n", nu, n, "GK"))
    } else {
      pr("GK", fg$par, fg$loglik)
      cat(sprintf("%17s GK - exact: scale %+.6f  logLik %+.6f\n", "",
                  fg$par[["s"]] - fe$par[["s"]],
                  fg$loglik - fe$loglik))
    }
  }
}

cat("\n\nC. The worst corner from A3: n = 2, latent scale 0.25.\n\n")
for (nu in c(2.5, 3, 5)) {
  d <- sim_tre(G = 40, n = 2, s = 0.25, nu = nu, seed = 7001)
  fl <- fit_laplace(d, nu = nu)
  fe <- fit_exact(d, nu = nu, K = 101,
                  start = c(fl$par[1], fl$par[2], log(fl$par[3]),
                            log(fl$par[4])))
  fl$obj$fn(fl$opt$par)
  anchor <- fl$obj$env$parList(fl$opt$par, fl$obj$env$last.par)
  anchor$log_nu <- NULL
  fg <- tryCatch(fit_gk(d, nu, anchor[c("beta", "log_sigma", "log_s",
                                        "b")]), error = function(e) NULL)
  cat(sprintf("nu %-5s exact s %8.5f  Laplace s %8.5f  GK s %s\n", nu,
              fe$par[["s"]], fl$par[["s"]],
              if (is.null(fg)) "   FAILED" else
                sprintf("%8.5f", fg$par[["s"]])))
}

sink()
cat("done\n")
