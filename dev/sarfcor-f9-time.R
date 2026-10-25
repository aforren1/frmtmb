# Fix round 1, F9: pearson residuals and autocor_matrix() of a sar() fit
# at n = 2500 (a 50 x 50 rook lattice, row-standardized), before and
# after the change. "before" puts the round's earlier code back in the
# loaded namespaces: autocor_seam_B() solving one column at a time,
# sar()'s unwhiten() on a vector only, and the pearson marginal read off
# the diagonal of B B'. Rscript <this> [before|after]
mode <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
k <- 50
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
nb <- lapply(seq_len(n), function(i) {
  which(abs(g$r - g$r[i]) + abs(g$c - g$c[i]) == 1)
})
W <- Matrix::sparseMatrix(i = rep(seq_len(n), lengths(nb)), j = unlist(nb),
                          x = rep(1 / lengths(nb), lengths(nb)),
                          dims = c(n, n))
set.seed(1)
d <- data.frame(x = rnorm(n))
A <- Matrix::Diagonal(n) - 0.5 * W
d$y <- as.vector(Matrix::solve(A, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
t_fit <- system.time(fit <- frm(y ~ x + sar(W), data = d,
                                data2 = list(W = W)))[["elapsed"]]
if (identical(mode, "before")) {
  ns <- asNamespace("frmtmb")
  old_B <- function(ac, theta, mu, sigma) {
    e <- ac_block_entry(ac)
    n <- ac[["d"]]
    mu <- rep_len(as.numeric(mu), n)
    sigma <- as.numeric(sigma)
    m <- e[["unwhiten"]](numeric(n), mu, sigma, theta, ac[["aux"]])
    vapply(seq_len(n), function(j) {
      u <- numeric(n)
      u[j] <- 1
      as.numeric(e[["unwhiten"]](u, mu, sigma, theta, ac[["aux"]])) - m
    }, numeric(n))
  }
  environment(old_B) <- ns
  old_marg <- function(ac, theta) diag(autocor_cor(theta, ac))
  environment(old_marg) <- ns
  utils::assignInNamespace("autocor_seam_B", old_B, ns = "frmtmb")
  utils::assignInNamespace("autocor_seam_marginal", old_marg, ns = "frmtmb")
  old_unwhiten <- function(u, mu, sigma, theta, aux) {
    rho <- aux[["lo"]] + (aux[["hi"]] - aux[["lo"]]) / (1 + exp(-theta[1]))
    A <- Matrix::Diagonal(nrow(aux[["W"]])) - rho * aux[["W"]]
    if (identical(aux[["type"]], "lag")) {
      as.vector(Matrix::solve(A, mu + sigma * u))
    } else {
      mu + as.vector(Matrix::solve(A, sigma * u))
    }
  }
  reg <- get("frmtmb_autocor_registry", ns)
  reg$reg$sar$unwhiten <- old_unwhiten
}
t_p <- system.time(rp <- residuals(fit, type = "pearson"))[["elapsed"]]
t_m <- system.time(M <- autocor_matrix(fit))[["elapsed"]]
saveRDS(list(rp = rp[, "Estimate"], diagM = diag(M), M12 = M[1:5, 1:5]),
        file.path(tempdir(), "x.rds"))
cat(sprintf("%s n %d | frm %.2f s | pearson %.2f s | autocor_matrix %.2f s | sum(rp) %.12g | sum(diag M) %.12g\n",
            mode, n, t_fit, t_p, t_m, sum(rp[, "Estimate"]), sum(diag(M))))
