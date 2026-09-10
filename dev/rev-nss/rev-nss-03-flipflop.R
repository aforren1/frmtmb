# REVIEW of lane nss, attack 2 continued, and attack 8.
#
# rev-nss-02-geom.R found 9 rows of 598 where the shipped default is
# WORSE than truncation, worst 2.10x. Every one is flip-flop: the
# absorption mode `exp(-ka * ii)` is close to the slow disposition mode
# `exp(-lambda_z * ii)`, and in the CENTRAL compartment those two modes
# carry opposite signs, which is the one arrangement that makes the
# measured ratio overstate the dominant mode instead of understating it.
#
# The lane's own 22-schedule sweep (nss-12-warn.R) cannot see this: every
# oral row there has `ka` at 1.0 or 1.1 per hour against a `lambda_z` of
# at most 0.2, so the absorption mode is 5x to 1000x faster than the
# terminal one and has decayed to nothing by cycle 18. The regime is not
# rare in the field; it is what an extended-release or depot formulation
# is written to produce.
#
# Two things have to be settled before this is a blocker rather than a
# curiosity:
#   A. WHERE IS THE BOUNDARY. If extrapolation only loses where
#      truncation was already useless, the finding is much weaker. So:
#      the smallest truncated error on any row where extrapolation
#      loses, against the largest truncated error on any row where it
#      wins.
#   B. DOES THE WARNING SAY SO. Attack 8's second case: the guarded
#      thing is PRESENT but the guard is quiet.
# And the instrument is checked in C: the worst row is run through
# `frm_ode()` itself, not through this file's cycle map.
#
# Script path: dev/rev-nss/rev-nss-03-flipflop.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()
ext <- frmtmb.ode:::ode_ss_extrapolate
SS_TOL <- 1e-6

build <- function(ke, ka, k12, k21, ii, dose = 100) {
  M <- matrix(c(-ka, 0, 0,
                ka, -(ke + k12), k21,
                0, k12, -k21), 3, 3, byrow = TRUE)
  # Matrix::expm, not an eigendecomposition: ka = lambda_z exactly is
  # the middle of the region under study and makes the eigenvectors
  # singular.
  A <- as.matrix(Matrix::expm(M * ii))
  e <- c(dose, 0, 0)
  list(A = A, e = e,
       yinf = as.numeric(solve(diag(3) - A, A %*% e)))
}
runin <- function(m, n_ss) {
  y <- numeric(3)
  keep <- vector("list", n_ss + 1L)
  keep[[1L]] <- y
  for (k in seq_len(n_ss)) {
    y <- as.numeric(m$A %*% (y + m$e))
    keep[[k + 1L]] <- y
  }
  keep[seq(n_ss - 1L, n_ss + 1L)]
}
one <- function(ke, ka, k12, k21, ii, n_ss = 20L,
                atol = 1e-8, rtol = 1e-8) {
  m <- build(ke, ka, k12, k21, ii)
  kp <- runin(m, n_ss)
  y <- kp[[3L]]
  ex <- ext(kp[[1L]], kp[[2L]], y, atol, rtol)
  sc <- max(1e-12, max(abs(m$yinf)))
  list(trunc = max(abs(y - m$yinf)) / sc,
       extrap = max(abs(as.numeric(ex$y) - m$yinf)) / sc,
       r = as.numeric(ex$r),
       # what ode_run_in_finish() reports on the default arm: the
       # disagreement between successive extrapolants
       rel = local({
         kp4 <- runin(m, n_ss)
         m2 <- build(ke, ka, k12, k21, ii)
         full <- local({
           yy <- numeric(3); kk <- vector("list", n_ss + 1L)
           kk[[1L]] <- yy
           for (k in seq_len(n_ss)) {
             yy <- as.numeric(m2$A %*% (yy + m2$e)); kk[[k + 1L]] <- yy
           }
           kk
         })
         a <- ext(full[[n_ss - 2L]], full[[n_ss - 1L]], full[[n_ss]],
                  atol, rtol)
         b <- ext(full[[n_ss - 1L]], full[[n_ss]], full[[n_ss + 1L]],
                  atol, rtol)
         max(abs(as.numeric(b$y) - as.numeric(a$y))) /
           max(1e-12, max(abs(as.numeric(b$y))))
       }))
}

lamz <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b - sqrt(b * b - 4 * ke * k21)) / 2
}

## ---- A. the boundary ------------------------------------------------
cat("\n== A. boundary of the geometric assumption, 2 cmt oral,",
    "n_ss = 20 ==\n")
rows <- list()
for (lzii in c(0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3)) {
  for (krat in c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2,
                 1.5, 2, 3, 5, 10, 50)) {
    # the peripheral shape is swept too: it sets the SIZE of the second
    # disposition mode, which is what decides whether the two dominant
    # weights can cancel
    for (sh in list(c(4, 8, 1.2), c(2, 2, 0.33), c(10, 30, 1.05),
                    c(1.5, 0.6, 1.5))) {
      ii <- 24
      lz <- lzii / ii
      ke <- sh[1L] * lz; k12 <- sh[2L] * lz; k21 <- sh[3L] * lz
      lzz <- lamz(ke, k12, k21)
      ke <- ke * lz / lzz; k12 <- k12 * lz / lzz; k21 <- k21 * lz / lzz
      ka <- krat * lz
      o <- one(ke, ka, k12, k21, ii)
      rows[[length(rows) + 1L]] <- data.frame(
        lzii = lzii, krat = krat, shape = paste(sh, collapse = "/"),
        trunc = o$trunc, extrap = o$extrap,
        ratio = o$extrap / o$trunc, rel = o$rel,
        warns = o$rel > SS_TOL)
    }
  }
}
d <- do.call(rbind, rows)
d <- d[is.finite(d$ratio), ]
lose <- d[d$ratio > 1.0000001, ]
win <- d[d$ratio <= 1.0000001, ]
cat("rows:", nrow(d), "  extrapolation loses on:", nrow(lose), "\n")
cat(sprintf("smallest truncated error on a LOSING row: %.3e\n",
            min(lose$trunc)))
cat(sprintf("largest truncated error on a WINNING row:  %.3e\n",
            max(win$trunc)))
cat(sprintf("worst ratio: %.2f (trunc %.3e -> extrap %.3e)\n",
            max(lose$ratio), lose$trunc[which.max(lose$ratio)],
            lose$extrap[which.max(lose$ratio)]))
cat("\nlosing rows, worst 12 by ratio:\n")
print(format(head(lose[order(-lose$ratio), ], 12), digits = 3),
      row.names = FALSE)
cat("\nka/lambda_z of every losing row:\n")
print(sort(unique(lose$krat)))

## ---- B. is the warning quiet on those rows? ------------------------
cat("\n== B. attack 8: the guarded thing PRESENT, is the guard quiet? ==\n")
cat("A row is a MISS if extrapolated error > ss_tol =", SS_TOL,
    "and rel <= ss_tol.\n\n")
d$miss <- d$extrap > SS_TOL & !d$warns
d$false_alarm <- d$extrap <= SS_TOL & d$warns
cat("misses:", sum(d$miss), " false alarms:", sum(d$false_alarm),
    " of", nrow(d), "rows\n")
if (any(d$miss))
  print(format(head(d[d$miss, ][order(-d[d$miss, ]$extrap), ], 12),
               digits = 3), row.names = FALSE)
cat("\nunderstatement rel/extrap on the LOSING rows:\n")
print(summary(lose$rel / lose$extrap))

## ---- C. the instrument: the worst row through frm_ode() ------------
cat("\n== C. same design through frm_ode(), not through this file ==\n")
lz <- 0.05 / 24; ke0 <- 4 * lz; k12 <- 8 * lz; k21 <- 1.2 * lz
s <- lz / lamz(ke0, k12, k21)
ke <- ke0 * s; k12 <- k12 * s; k21 <- k21 * s
ka <- 1.5 * lz
cat(sprintf("ke %.6f  k12 %.6f  k21 %.6f  ka %.6f  lambda_z %.6f",
            ke, k12, k21, ka, lamz(ke, k12, k21)),
    sprintf(" t1/2z %.1f h  ii 24\n", log(2) / lamz(ke, k12, k21)))
dyn <- function(t, y, p) list(c(-p[4L] * y[1L],
                                p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] +
                                  p[3L] * y[3L],
                                p[2L] * y[2L] - p[3L] * y[3L]))
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
  list(v = as.numeric(v) / 10, w = w)
}
exact <- as.numeric(frm_lincmt(
  parms = list(ke = ke, k12 = k12, k21 = k21, ka = ka, V = 10),
  times = tt, ncmt = 2, depot = TRUE,
  events = data.frame(time = 0, state = "depot", value = 100, ii = 24,
                      addl = 0L, ss = TRUE)))
tr <- run_ode(FALSE); ex2 <- run_ode(TRUE)
sc <- max(abs(exact))
cat(sprintf("truncated   max rel err over the interval: %.3e\n",
            max(abs(tr$v - exact)) / sc))
cat(sprintf("extrapolated max rel err over the interval: %.3e\n",
            max(abs(ex2$v - exact)) / sc))
cat(sprintf("ratio: %.2f\n",
            max(abs(ex2$v - exact)) / max(abs(tr$v - exact))))
cat("\nwarning, ss_extrapolate = FALSE:\n  ",
    if (is.null(tr$w)) "(none)" else paste(tr$w, collapse = "\n  "), "\n")
cat("\nwarning, ss_extrapolate = TRUE (the shipped default):\n  ",
    if (is.null(ex2$w)) "(none)" else paste(ex2$w, collapse = "\n  "),
    "\n")
cat("\ndone\n")
