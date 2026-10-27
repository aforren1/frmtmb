# Reviewer of lane loglik: the whole per-draw path of frmtmb.sample's
# log_lik() (draws_fit_at(), eval_dpars(), the dispatch) against
# autocor_row_loglik() alone, 100 draws, min of 3 rounds. Seed 95.
# Output: dev/review-loglik-log/e2e.txt
LIB <- Sys.getenv("RV_LIB", "C:/Users/adf44/source/r/wt-loglik-lib")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
draws_at <- function(fit, vals) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  k <- max(lengths(vals))
  m <- matrix(0, k, length(lab), dimnames = list(NULL, lab))
  for (nm in names(vals)) m[, nm] <- vals[[nm]]
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}
set.seed(95)
K <- 100L
timeit <- function(f) min(vapply(1:3, function(i) {
  system.time(f())[["elapsed"]]
}, 0))

# sar() on a 40 x 40 rook lattice
k <- 40
n <- k * k
id <- matrix(seq_len(n), k, k)
P <- matrix(0, n, n)
P[cbind(as.vector(id[-k, ]), as.vector(id[-1, ]))] <- 1
P[cbind(as.vector(id[, -k]), as.vector(id[, -1]))] <- 1
P <- P + t(P)
W <- P / rowSums(P)
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n)
fit <- frm(y ~ x + sar(W), data = d, data2 = list(W = W),
           dry_run = "objective")
aux <- fit$frame$autocor$y$aux
rho <- runif(K, 0.1, 0.6)
ds <- draws_at(fit, list(b_Intercept = rnorm(K, 1, 0.1),
                         b_x = rnorm(K, 0.5, 0.1),
                         sigma = exp(rnorm(K, 0, 0.05)),
                         thetaac_1 = qlogis((rho - aux$lo) /
                                              (aux$hi - aux$lo))))
idx <- frmtmb.sample:::draws_index(ds)
fits <- lapply(seq_len(K), function(s) frmtmb.sample:::draws_fit_at(ds, s, idx))
dps <- lapply(fits, eval_dpars)
t_all <- timeit(function() log_lik(ds))
t_core <- timeit(function() for (s in seq_len(K)) {
  autocor_row_loglik(fits[[s]], dps[[s]], "y")
})
cat(sprintf("sar n = %d, %d draws: log_lik() %.3f s, autocor_row_loglik() alone %.3f s, overhead share %.2f\n",
            n, K, t_all, t_core, 1 - t_core / t_all))

# ar(cov = TRUE), 400 groups x 10
G <- 400
d2 <- data.frame(g = rep(seq_len(G), each = 10), t = rep(1:10, G))
d2$x <- rnorm(nrow(d2))
d2$y <- rnorm(nrow(d2))
fit2 <- frm(y ~ x + ar(t, g, cov = TRUE), data = d2, dry_run = "objective")
a <- runif(K, 0.2, 0.7)
ds2 <- draws_at(fit2, list(b_Intercept = rnorm(K), b_x = rnorm(K),
                           sigma = exp(rnorm(K, 0, 0.05)),
                           thetaac_1 = a / sqrt(1 - a^2)))
idx2 <- frmtmb.sample:::draws_index(ds2)
fits2 <- lapply(seq_len(K), function(s) {
  frmtmb.sample:::draws_fit_at(ds2, s, idx2)
})
dps2 <- lapply(fits2, eval_dpars)
t_all2 <- timeit(function() log_lik(ds2))
t_core2 <- timeit(function() for (s in seq_len(K)) {
  autocor_row_loglik(fits2[[s]], dps2[[s]], "y")
})
cat(sprintf("ar 400 x 10, %d draws: log_lik() %.3f s, autocor_row_loglik() alone %.3f s, overhead share %.2f\n",
            K, t_all2, t_core2, 1 - t_core2 / t_all2))
cat("DONE\n")
