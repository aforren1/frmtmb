# Lane predraw2: which proposal should the reweighting draw from where the
# check fires? On the sweep's three weakly identified mixtures (seed 101,
# or the seeds given), importance sampling likelihood / proposal with
#   norm    N(theta_hat, V), the 0.82.0 draw law (N proposals)
#   t4      multivariate t, 4 df, scale V
#   adapt   one adaptation: normal with the PSIS-weighted mean and
#           covariance of the norm stage, then t4 around those
# each resampled to 1000 parameter draws, and the predictive simulated at
# each draw. Reported: Pareto k, ESS share, unique draws kept, and the
# predictive's coverage of the true law (dev/predraw2-mixcov.R's truth).
# Usage: Rscript dev/predraw2-prop.R <lib> <design> <seeds> <N>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
design <- strsplit(args[2], ":")[[1]]
seeds <- eval(parse(text = args[3]))
N <- as.integer(args[4])
fam <- design[1]
case <- design[2]

objdrop <- function(fit, D) {
  obj <- fit$obj
  p <- fit$opt$par
  saved <- ns$obj_state_save(obj)
  on.exit(ns$obj_state_restore(obj, saved))
  f0 <- obj$fn(p)
  off <- nrow(D) - length(p)
  vapply(seq_len(ncol(D)), function(s) {
    ns$obj_state_restore(obj, saved)
    f1 <- tryCatch(suppressWarnings(obj$fn(D[off + seq_along(p), s])),
                   error = function(e) NA_real_)
    if (is.finite(f1)) f1 - f0 else Inf
  }, 0)
}
rmvt <- function(N, mu, S, df) {
  L <- chol(S)
  Z <- matrix(stats::rnorm(length(mu) * N), length(mu), N)
  w <- sqrt(df / stats::rchisq(N, df))
  X <- crossprod(L, Z) * rep(w, each = length(mu))
  list(D = mu + X, lq = -(df + length(mu)) / 2 *
         log1p(colSums(backsolve(L, X, transpose = TRUE)^2) / df))
}
rmvn <- function(N, mu, S) {
  L <- chol(S)
  Z <- matrix(stats::rnorm(length(mu) * N), length(mu), N)
  X <- crossprod(L, Z)
  list(D = mu + X, lq = -colSums(Z^2) / 2)
}
isr <- function(fit, pr) {
  lw <- -objdrop(fit, pr$D) - pr$lq
  ps <- suppressWarnings(loo::psis(lw, r_eff = 1))
  w <- exp(as.numeric(ps$log_weights) - max(ps$log_weights))
  w <- w / sum(w)
  list(k = ps$diagnostics$pareto_k, ess = 1 / sum(w^2) / length(w), w = w,
       D = pr$D)
}
resample <- function(w, M) {
  cw <- cumsum(w)
  cw[length(cw)] <- 1
  findInterval((stats::runif(1) + seq_len(M) - 1) / M, cw) + 1L
}
sim_at <- function(fit, D, map) {
  t(vapply(seq_len(ncol(D)), function(s) {
    fs <- ns$fit_set_outer(fit, D[, s], map)
    as.numeric(suppressWarnings(predict(fs, ndraws = 1L, summary = FALSE,
                                        propagate_error = FALSE)))
  }, numeric(n)))
}
src <- readLines("dev/predraw2-mixcov.R")
eval(parse(text = src[grep("^truth_row <-", src):
                        (grep("^rows <- list", src) - 1L)]))
for (seed in seeds) {
  d <- gen(fam, case, seed)
  s <- spec(fam)
  fit <- suppressMessages(suppressWarnings(frm(s[[1]], family = s[[2]],
                                               data = d)))
  ds <- ns$fit_draw_space(fit)
  V <- ds$V
  map <- ds$map
  v0 <- ns$fit_outer_vector(fit, map)
  tr <- truth_row(d)
  set.seed(1)
  arms <- list()
  a <- isr(fit, rmvn(N, v0, V))
  arms$norm <- a
  arms$t4 <- isr(fit, rmvt(N, v0, V, 4))
  mu <- colSums(t(a$D) * a$w)
  Xc <- a$D - mu
  S <- (Xc * rep(a$w, each = nrow(Xc))) %*% t(Xc)
  S <- (S + t(S)) / 2 + diag(1e-8 * diag(V), nrow(V))
  arms$adapt <- tryCatch(isr(fit, rmvt(N, mu, S, 4)),
                         error = function(e) NULL)
  for (nm in names(arms)) {
    a <- arms[[nm]]
    if (is.null(a)) next
    idx <- resample(a$w, 1000L)
    p <- sim_at(fit, a$D[, idx, drop = FALSE], map)
    if (tr$type == "cat") {
      sh <- t(apply(p, 2L, function(col) tabulate(col, 4L) / 1000))
      met <- sprintf("mae %.4f", mean(abs(sh - tr$P)))
    } else {
      q <- apply(p, 2L, stats::quantile, c(0.01, 0.05, 0.95, 0.99),
                 names = FALSE)
      c90 <- mean(vapply(seq_len(n), function(i) mass_in(tr, i, q[2, i],
                                                         q[3, i]), 0))
      c98 <- mean(vapply(seq_len(n), function(i) mass_in(tr, i, q[1, i],
                                                         q[4, i]), 0))
      met <- sprintf("cov90 %.3f cov98 %.3f w98 %.3f", c90, c98,
                     mean(q[4, ] - q[1, ]))
    }
    cat(sprintf("seed %3d %-6s N %5d k %5.2f ess %.4f unique %4d %s\n", seed,
                nm, N, a$k, a$ess, length(unique(idx)), met))
  }
}
