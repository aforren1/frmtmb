# Lane loglik: the cost of one draw's pointwise log-density for a
# residual covariance term, measured. Arms, interleaved in one process,
# each block repeated past 1.2 s, minimum of 5 rounds:
#   sar   : frmtmb's route (sparse Jacobian (I - rho W) / sigma)
#   sarfb : the same term with no jacobian slot (the fallback, a dense
#           inverse of unwhiten() of the identity), registered under
#           another name
#   brms  : brms 2.23.0's algorithm for the same draw
#           (log_lik_gaussian_lagsar(): dense I - rho W, its crossproduct
#           and a dense solve), transcribed
#   ctrl  : the sar arm again, which must report 1.0
# plus ar(cov = TRUE) on 400 groups of 10 rows (one pattern) and fcor()
# on 600 rows. Allocation: bytes from Rprofmem() per call when R has
# memory profiling, otherwise the change in gc()'s max-used cells.
# Seed 5. Output: dev/loglik-cost.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spatial)})
cat("lib:", find.package("frmtmb"), "\n")
cat("profmem:", capabilities("profmem"), "\n")
set.seed(5)
k <- 40
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * d$x + rnorm(n)))
fit <- frm(y ~ x + sar(W), data = d, data2 = list(W = W))
cat("sar n =", n, " nnz(W) =", sum(W != 0), "\n")

# the fallback arm: sar's own map without its jacobian slot
frmtmb_local_registries()
e <- frmtmb:::ac_entry("sar")
frmtmb_register_autocor("sarfb", args = c("M", "type"), parse = e$parse,
                        build = e$build, whiten = e$whiten,
                        unwhiten = e$unwhiten, npar = 1L, start = e$start,
                        natural = e$natural, sigma_varies = FALSE,
                        bounds = e$bounds)
ffb <- frm(y ~ x + sarfb(W), data = d, data2 = list(W = W))
dpv <- eval_dpars(fit)
dpb <- eval_dpars(ffb)
rho <- frmtmb:::autocor_natural(fit$estimates$thetaac, fit$frame$autocor$y)
mu <- as.numeric(frm_linpred(fit))
sig <- sigma(fit)
brms_lagsar <- function() {
  IB <- diag(n) - rho * W
  Cinv <- t(IB) %*% IB / sig^2
  e <- d$y - solve(IB, mu)
  gg <- Cinv %*% e
  cbar <- diag(Cinv)
  dnorm(d$y, d$y - gg / cbar, sqrt(1 / cbar), log = TRUE)
}
arms <- list(
  sar = function() autocor_row_loglik(fit, dpv, "y"),
  sarfb = function() autocor_row_loglik(ffb, dpb, "y"),
  brms = brms_lagsar,
  ctrl = function() autocor_row_loglik(fit, dpv, "y"))
v <- lapply(arms, function(f) as.numeric(f()))
cat(sprintf("agreement: sar vs brms %.3e, sarfb vs sar %.3e (max abs / max)\n",
            max(abs(v$sar - v$brms)) / max(abs(v$brms)),
            max(abs(v$sarfb - v$sar)) / max(abs(v$sar))))

time_arms <- function(arms, rounds = 5) {
  reps <- vapply(arms, function(f) {
    r <- 1L
    repeat {
      t0 <- proc.time()[["elapsed"]]
      for (i in seq_len(r)) f()
      el <- proc.time()[["elapsed"]] - t0
      if (el > 1.2) break
      r <- r * 2L
    }
    r
  }, 1L)
  best <- stats::setNames(rep(Inf, length(arms)), names(arms))
  for (k in seq_len(rounds)) {
    for (nm in names(arms)) {
      t0 <- proc.time()[["elapsed"]]
      for (i in seq_len(reps[[nm]])) arms[[nm]]()
      best[[nm]] <- min(best[[nm]],
                        (proc.time()[["elapsed"]] - t0) / reps[[nm]])
    }
  }
  list(best = best, reps = reps)
}
alloc <- function(f) {
  if (capabilities("profmem")) {
    tf <- tempfile()
    utils::Rprofmem(tf, threshold = 0)
    f()
    utils::Rprofmem(NULL)
    ln <- readLines(tf)
    b <- suppressWarnings(as.numeric(sub(":.*", "", ln)))
    return(sum(b, na.rm = TRUE))
  }
  NA_real_
}
report <- function(tm, arms) {
  for (nm in names(arms)) {
    cat(sprintf("  %-6s %9.3f ms per draw (reps %4d)  alloc %10.0f bytes\n",
                nm, 1e3 * tm$best[[nm]], tm$reps[[nm]], alloc(arms[[nm]])))
  }
}
tm <- time_arms(arms)
cat("sar, per draw, min over 5 rounds:\n")
report(tm, arms)
cat(sprintf("  brms / sar = %.1f, sarfb / sar = %.1f, ctrl / sar = %.2f\n",
            tm$best[["brms"]] / tm$best[["sar"]],
            tm$best[["sarfb"]] / tm$best[["sar"]],
            tm$best[["ctrl"]] / tm$best[["sar"]]))

# ar(cov = TRUE): 400 groups of 10
G <- 400
da <- expand.grid(t = 1:10, gr = factor(seq_len(G)))
da$x <- rnorm(nrow(da))
da$y <- 1 + 0.5 * da$x + as.vector(replicate(G, as.numeric(
  stats::arima.sim(list(ar = 0.6), 10))))
fa <- frm(y ~ x + ar(t, gr, cov = TRUE), data = da)
dpa <- eval_dpars(fa)
R <- autocor_matrix(fa)
sa <- sigma(fa)
mua <- as.numeric(frm_linpred(fa))
brms_time <- function() {
  # brms's log_lik_gaussian_time(): per group, a dense solve and an
  # inverse of the group's covariance
  out <- numeric(nrow(da))
  for (gg in levels(da$gr)) {
    obs <- which(da$gr == gg)
    C <- sa^2 * R
    Cinv <- solve(C)
    e <- da$y[obs] - mua[obs]
    gv <- solve(C, e)
    cbar <- diag(Cinv)
    out[obs] <- dnorm(da$y[obs], da$y[obs] - gv / cbar, sqrt(1 / cbar),
                      log = TRUE)
  }
  out
}
arms2 <- list(ar = function() autocor_row_loglik(fa, dpa, "y"),
              brms = brms_time,
              ctrl = function() autocor_row_loglik(fa, dpa, "y"))
cat(sprintf("ar agreement with brms's algorithm: %.3e\n",
            max(abs(arms2$ar() - brms_time())) / max(abs(brms_time()))))
tm <- time_arms(arms2)
cat("ar(cov = TRUE), 400 groups x 10, per draw:\n")
report(tm, arms2)
cat(sprintf("  brms / ar = %.1f, ctrl / ar = %.2f\n",
            tm$best[["brms"]] / tm$best[["ar"]],
            tm$best[["ctrl"]] / tm$best[["ar"]]))

# fcor(): 600 rows
nf <- 600
Vf <- 0.7^abs(outer(1:nf, 1:nf, "-"))
df <- data.frame(x = rnorm(nf))
df$y <- 1 + 0.5 * df$x + as.vector(t(chol(Vf)) %*% rnorm(nf))
ff <- frm(y ~ x + fcor(Vf), data = df, data2 = list(Vf = Vf))
dpf <- eval_dpars(ff)
sf <- sigma(ff)
muf <- as.numeric(frm_linpred(ff))
brms_fcor <- function() {
  C <- sf^2 * Vf
  Cinv <- solve(C)
  e <- df$y - muf
  gv <- solve(C, e)
  cbar <- diag(Cinv)
  dnorm(df$y, df$y - gv / cbar, sqrt(1 / cbar), log = TRUE)
}
arms3 <- list(fcor = function() autocor_row_loglik(ff, dpf, "y"),
              brms = brms_fcor,
              ctrl = function() autocor_row_loglik(ff, dpf, "y"))
cat(sprintf("fcor agreement with brms's algorithm: %.3e\n",
            max(abs(arms3$fcor() - brms_fcor())) / max(abs(brms_fcor()))))
tm <- time_arms(arms3)
cat("fcor(), 600 rows, per draw:\n")
report(tm, arms3)
cat(sprintf("  brms / fcor = %.1f, ctrl / fcor = %.2f\n",
            tm$best[["brms"]] / tm$best[["fcor"]],
            tm$best[["ctrl"]] / tm$best[["fcor"]]))
