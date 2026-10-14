# The simplex families' log densities against brms 2.23.0's own Stan
# functions (inst/chunks of the brms tarball in dev/brms-suite), exposed
# through rstan::expose_stan_functions(), on a grid.
# Output: dev/famsimplex-stan-lpdf.log. Seed 20261013.
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(rstan)})
cat("lib:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex"
chunks <- file.path(wt, "dev/brms-suite/brms/inst/chunks")
code <- paste0(
  "functions {\n",
  paste(readLines(file.path(chunks, "fun_dirichlet_logit.stan")),
        collapse = "\n"), "\n",
  paste(readLines(file.path(chunks, "fun_logistic_normal.stan")),
        collapse = "\n"), "\n",
  # wrappers callable from R; the dirichlet2 likelihood is Stan's own
  # dirichlet_lpdf, which brms's stan_log_lik_dirichlet2() writes
  "  real w_dirichlet_logit(vector y, vector mu, real phi) {\n",
  "    return dirichlet_logit_lpdf(y | mu, phi);\n  }\n",
  "  real w_dirichlet2(vector y, vector alpha) {\n",
  "    return dirichlet_lpdf(y | alpha);\n  }\n",
  "  real w_logistic_normal(vector y, vector mu, vector sigma, matrix Lcor,\n",
  "                         int ref) {\n",
  "    return logistic_normal_cholesky_cor_lpdf(y | mu, sigma, Lcor, ref);\n",
  "  }\n}\n")
env <- new.env()
expose_stan_functions(stanc(model_code = code), env = env)
cat("stan functions:", paste(ls(env), collapse = " "), "\n")

ns <- asNamespace("frmtmb")
set.seed(20261013)
rdir <- function(n, a) {
  G <- matrix(rgamma(n * length(a), a), n, byrow = TRUE)
  G / rowSums(G)
}
rcorr <- function(D) {
  W <- crossprod(matrix(rnorm(D * (D + 2)), D + 2, D))
  stats::cov2cor(W)
}
res <- list()
add <- function(fam, K, ref, a, b) {
  d <- abs(a - b)
  res[[length(res) + 1L]] <<- data.frame(
    family = fam, K = K, ref = ref, n = length(a),
    max_abs = max(d), max_rel = max(d / pmax(abs(b), 1)),
    all_finite = all(is.finite(a)) && all(is.finite(b)))
}
for (K in 2:5) {
  cats <- letters[seq_len(K)]
  for (ref in c(seq_len(K), NA)) {
    n <- 200
    # concentrations from tiny to large, so components near 0 occur
    Y <- rdir(n, rep(c(0.3, 1, 5)[(K %% 3) + 1L], K))
    Y <- pmax(Y, 1e-12)
    Y <- Y / rowSums(Y)
    # the dirichlet family, reference at ref (NA: every category)
    fam <- ns$fam_dirichlet_impl(cats, if (is.na(ref)) NA_integer_ else ref)
    sx <- fam$simplex
    E <- matrix(rnorm(n * length(sx$pred), 0, 2), n)
    E[1:5, ] <- 30 * sign(E[1:5, ])
    phi <- exp(runif(n, log(0.05), log(1e4)))
    dp <- c(stats::setNames(lapply(seq_along(sx$pred), function(j) E[, j]),
                            sx$mu), list(phi = phi))
    ours <- fam$lpdf(Y, dp, list())
    stan <- vapply(seq_len(n), function(i) {
      mu <- numeric(K)
      mu[sx$pred] <- E[i, ]
      env$w_dirichlet_logit(Y[i, ], mu, phi[i])
    }, 0)
    add("dirichlet", K, ref, ours, stan)
    if (is.na(ref)) {
      # dirichlet2: every category its own concentration
      f2 <- ns$fam_dirichlet2_impl(cats)
      A <- matrix(exp(runif(n * K, log(0.05), log(500))), n)
      dp2 <- stats::setNames(lapply(seq_len(K), function(k) A[, k]),
                             f2$simplex$mu)
      ours2 <- f2$lpdf(Y, dp2, list())
      stan2 <- vapply(seq_len(n), function(i) env$w_dirichlet2(Y[i, ], A[i, ]),
                      0)
      add("dirichlet2", K, NA, ours2, stan2)
      next
    }
    # the logistic-normal family
    fl <- ns$fam_logistic_normal_impl(cats, ref)
    sl <- fl$simplex
    D <- K - 1L
    M <- matrix(rnorm(n * D, 0, 2), n)
    S <- matrix(exp(runif(n * D, log(0.05), log(5))), n)
    C <- if (D > 1L) rcorr(D) else matrix(1, 1, 1)
    dpl <- c(stats::setNames(lapply(seq_len(D), function(j) M[, j]), sl$mu),
             stats::setNames(lapply(seq_len(D), function(j) S[, j]),
                             sl$sigma))
    ex <- if (D > 1L) list(lncor_raw = ns$us_theta_cor(C)) else list()
    oursl <- fl$lpdf(Y, dpl, list(), ex)
    Lc <- t(chol(C))
    stanl <- vapply(seq_len(n), function(i) {
      env$w_logistic_normal(Y[i, ], M[i, ], S[i, ], Lc, ref)
    }, 0)
    add("logistic_normal", K, ref, oursl, stanl)
  }
}
out <- do.call(rbind, res)
print(out, digits = 3, row.names = FALSE)
cat("\nWORST per family:\n")
print(aggregate(cbind(max_abs, max_rel) ~ family, out, max), digits = 3)
cat("rows compared:", sum(out$n), " any non-finite:",
    any(!out$all_finite), "\n")
