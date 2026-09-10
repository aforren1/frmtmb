# NIT 3: extrapolation loses to truncation on flip-flop absorption,
# which is ordinary linear pharmacokinetics and not a constructed
# oscillator.
#
# The mechanism is the reviewer's: the ratio the run-in reads is a
# weighted mean of the cycle map's modes, and when two weights have
# opposite signs the mean lies OUTSIDE their range and the correction
# overshoots. Opposite signs are the normal arrangement in an oral
# model, so the hazard is `ka` near `lambda_z`, which is what an
# extended-release or depot formulation is built to produce.
#
# Section A is the reviewer's `dev/rev-nss/rev-nss-03-flipflop.R`
# section C reproduced exactly, unrounded, against the SHIPPED shape.
# Section B sweeps for the bound: how bad can the loss get, and how
# wrong is truncation wherever it happens.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-23-flipflop.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)
                  library(Matrix)})
nss_report_env()

lamz <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b - sqrt(b * b - 4 * ke * k21)) / 2
}
dyn <- function(t, y, p) list(c(-p[4L] * y[1L],
                                p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] +
                                  p[3L] * y[3L],
                                p[2L] * y[2L] - p[3L] * y[3L]))

cat("\n== A. the reviewer's worst row, through frm_ode() ==\n")
lz <- 0.05 / 24
ke0 <- 4 * lz; k12 <- 8 * lz; k21 <- 1.2 * lz
s <- lz / lamz(ke0, k12, k21)
ke <- ke0 * s; k12 <- k12 * s; k21 <- k21 * s
ka <- 1.5 * lz
cat(sprintf("ke %.6f  k12 %.6f  k21 %.6f  ka %.6f  lambda_z %.6f",
            ke, k12, k21, ka, lamz(ke, k12, k21)))
cat(sprintf("  t1/2z %.1f h  ii 24\n", log(2) / lamz(ke, k12, k21)))
tt <- seq(0, 24, length.out = 25)
ev <- data.frame(time = 0, state = 1L, value = 100, ii = 24, ss = TRUE)
run_ode <- function(x) {
  w <- NULL
  v <- withCallingHandlers(
    frm_ode(dyn, init = rep(list(0), 3L), times = tt,
            parms = list(ke, k12, k21, ka), events = ev, output = 2L,
            n_ss = 20L, ss_extrapolate = x, atol = 1e-10, rtol = 1e-10),
    warning = function(e) { w <<- c(w, conditionMessage(e))
                            invokeRestart("muffleWarning") })
  list(v = as.numeric(v) / 10, warned = any(grepl("run-in", w)))
}
exact <- as.numeric(frm_lincmt(
  parms = list(ke = ke, k12 = k12, k21 = k21, ka = ka, V = 10),
  times = tt, ncmt = 2, depot = TRUE,
  events = data.frame(time = 0, state = "depot", value = 100, ii = 24,
                      addl = 0L, ss = TRUE)))
tr <- run_ode(FALSE); ex <- run_ode(TRUE)
sc <- max(abs(exact))
et <- max(abs(tr$v - exact)) / sc
ee <- max(abs(ex$v - exact)) / sc
cat(sprintf("truncated    %.4e  (warned %s)\n", et, tr$warned))
cat(sprintf("extrapolated %.4e  (warned %s)\n", ee, ex$warned))
cat(sprintf("ratio %.2f\n", ee / et))

cat("\n== B. the bound, over exact cycle maps ==\n")
Amat <- function(ke, k12, k21, ka)
  matrix(c(-ka, 0, 0, ka, -(ke + k12), k21, 0, k12, -k21), 3L, 3L,
         byrow = TRUE)
# the SHIPPED function itself rather than a replica, so the sweep
# cannot drift away from what is installed
shipped <- function(ya, yb, yc, atol, rtol)
  frmtmb.ode:::ode_ss_extrapolate(ya, yb, yc, atol, rtol)[["y"]]
cyc <- function(ke, k12, k21, ka, ii, n) {
  M <- as.matrix(expm(Amat(ke, k12, k21, ka) * ii))
  a <- c(100, 0, 0); y <- c(0, 0, 0)
  ys <- matrix(0, n + 1L, 3L)
  for (k in seq_len(n)) { y <- as.numeric(M %*% (y + a))
                          ys[k + 1L, ] <- y }
  list(ys = ys, lim = as.numeric(solve(diag(3) - M, M %*% a)))
}
res <- NULL
for (lzii in c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3)) {
  for (kar in c(0.1, 0.3, 0.5, 0.8, 1, 1.5, 2, 3, 5, 10, 50)) {
    for (sh in list(c(4, 8, 1.2), c(2, 3, 0.5), c(10, 20, 3),
                    c(1.5, 1, 0.2))) {
      lzv <- lzii / 24
      k <- sh * lzv
      sc2 <- lzv / lamz(k[1], k[2], k[3])
      k <- k * sc2
      kav <- kar * lzv
      z <- cyc(k[1], k[2], k[3], kav, 24, 21L)
      e_t <- max(abs(z$ys[21L, ] - z$lim)) / max(abs(z$lim))
      ye <- shipped(z$ys[19L, ], z$ys[20L, ], z$ys[21L, ], 1e-8, 1e-8)
      e_e <- max(abs(ye - z$lim)) / max(abs(z$lim))
      res <- rbind(res, data.frame(lzii = lzii, kar = kar,
                                   trunc = e_t, ext = e_e,
                                   ratio = e_e / e_t))
    }
  }
}
bad <- !is.finite(res$ratio)
cat(sprintf("  %d cycle maps swept, %d dropped as non-finite\n",
            nrow(res), sum(bad)))
res <- res[!bad, ]
lose <- res[res$ratio > 1, ]
cat(sprintf("  extrapolation loses on %d of them, worst ratio %.2f\n",
            nrow(lose), if (nrow(lose)) max(lose$ratio) else 1))
if (nrow(lose)) {
  cat(sprintf("  smallest truncated error on a LOSING row: %.4e\n",
              min(lose$trunc)))
  cat(sprintf("  largest truncated error on a WINNING row: %.4e\n",
              max(res$trunc[res$ratio <= 1])))
  cat(sprintf("  lambda_z * ii on the losing rows: %s\n",
              paste(sort(unique(lose$lzii)), collapse = ", ")))
  cat(sprintf("  ka / lambda_z on the losing rows: %s\n",
              paste(sort(unique(lose$kar)), collapse = ", ")))
  cat("  worst rows:\n")
  print(head(lose[order(-lose$ratio), ], 6))
}
cat(sprintf("  median ratio over all rows: %.4f\n", median(res$ratio)))
