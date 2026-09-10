# REVIEW of lane nss, attacks 3 and 4.
#
# 3. "It costs no extra solve, 21 per group at n_ss = 20 either way."
#    Counted two ways: deSolve::lsoda entries (the lane's count, which
#    is load-independent and repeated here), and the AD TAPE, because
#    arithmetic on tape variables is not free even when solves are.
# 4. The gradient. `d/dlog(k21)` from the truncated tape at +10.37
#    where the exact limit is -1.30, the WRONG SIGN.
#
# Script path: dev/rev-nss/rev-nss-04-cost-grad.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
rev_env()

two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}
II <- 24
TT <- c(0, 2, 6, 12, 23.9)
EV <- data.frame(time = 0, state = 1L, value = 100, ii = II, ss = TRUE)
EVL <- data.frame(time = 0, state = "depot", value = 100, ii = II,
                  addl = 0L, ss = TRUE)

f_ode <- function(th, ext, n = 20L) {
  sum(frm_ode(two_oral, init = list(0, 0, 0), times = TT,
              parms = list(exp(th[1L]), exp(th[2L]), exp(th[3L]),
                           exp(th[4L])),
              events = EV, output = 2L, n_ss = n, ss_tol = 1,
              ss_extrapolate = ext, atol = 1e-10, rtol = 1e-10))
}
f_exact <- function(th)
  10 * sum(frm_lincmt(parms = list(ke = exp(th[1L]), k12 = exp(th[2L]),
                                   k21 = exp(th[3L]), ka = exp(th[4L]),
                                   V = 10),
                      times = TT, ncmt = 2, depot = TRUE, events = EVL))
fd <- function(f, x, h = 1e-5)
  vapply(seq_along(x), function(j) {
    a <- x; b <- x; a[j] <- a[j] + h; b[j] <- b[j] - h
    (f(a) - f(b)) / (2 * h)
  }, 0)

x0 <- log(c(0.15, 0.3, 0.02, 1.0))

## ---- 3a. the solve count, counted rather than asserted -------------
cat("\n== 3a. deSolve::lsoda entries per group ==\n")
count_solves <- function(ext, n) {
  # trace() with a tracer expression cannot reach a counter in this
  # frame, so the count goes through a binding replacement instead
  cnt <- 0L
  ns <- asNamespace("deSolve")
  orig <- get("lsoda", envir = ns)
  unlockBinding("lsoda", ns)
  assign("lsoda", function(...) { cnt <<- cnt + 1L; orig(...) },
         envir = ns)
  on.exit({ assign("lsoda", orig, envir = ns); lockBinding("lsoda", ns) },
          add = TRUE)
  invisible(f_ode(x0, ext, n))
  cnt
}
for (n in c(5L, 20L)) for (ext in c(FALSE, TRUE))
  cat(sprintf("  n_ss = %2d  ss_extrapolate = %-5s  lsoda entries: %d\n",
              n, ext, count_solves(ext, n)))

## ---- 3b. the tape, which the lane did not count --------------------
cat("\n== 3b. AD tape size, both arms, same design ==\n")
tape_of <- function(ext, n = 20L)
  MakeTape(function(th) {
    "c" <- RTMB::ADoverload("c")
    f_ode(th, ext, n)
  }, x0)
# A node count is load-independent, which a clock is not. The printed
# operation stack is one line per node, so its length IS the node count.
n_nodes <- function(tp)
  length(capture.output(tp$print(depth = 1))) - 1L
# What one call to the correction costs on its own, so the delta above
# has something to be compared against.
solo <- n_nodes(MakeTape(function(x) sum(x), rep(1, 9)))
withc <- n_nodes(MakeTape(function(x) {
  "c" <- RTMB::ADoverload("c")
  z <- frmtmb.ode:::ode_ss_extrapolate(x[1:3], x[4:6], x[7:9], 1e-8,
                                       1e-8)
  sum(z$y) + sum(z$r)
}, rep(1, 9)))
cat("  one ode_ss_extrapolate() on 3 states:", withc - solo,
    "nodes
")
for (n in c(4L, 5L, 20L, 40L)) {
  a <- n_nodes(tape_of(FALSE, n))
  b <- n_nodes(tape_of(TRUE, n))
  cat(sprintf("  n_ss = %2d  nodes truncated %6d  extrapolated %6d",
              n, a, b))
  cat(sprintf("  delta %+d (%+.3f%%)\n", b - a, 100 * (b - a) / a))
}

## ---- 4. the gradient ----------------------------------------------
cat("\n== 4. the gradient at the 107 h design, n_ss = 20 ==\n")
g_ex <- fd(f_exact, x0)
nm <- c("d/dlke", "d/dlk12", "d/dlk21", "d/dlka")
cat(sprintf("  %-14s %12s %12s %12s %12s\n", "arm", nm[1], nm[2],
            nm[3], nm[4]))
cat(sprintf("  %-14s %12.4f %12.4f %12.4f %12.4f\n", "EXACT limit (fd)",
            g_ex[1], g_ex[2], g_ex[3], g_ex[4]))
for (ext in c(FALSE, TRUE)) {
  tp <- tape_of(ext)
  g_ad <- as.numeric(tp$jacfun()(x0))
  g_fd <- fd(function(th) f_ode(th, ext), x0)
  cat(sprintf("  %-14s %12.4f %12.4f %12.4f %12.4f\n",
              paste0("AD ext=", ext), g_ad[1], g_ad[2], g_ad[3],
              g_ad[4]))
  cat(sprintf("  %-14s %12.4f %12.4f %12.4f %12.4f\n",
              paste0("fd ext=", ext), g_fd[1], g_fd[2], g_fd[3],
              g_fd[4]))
  cat(sprintf("     AD vs its own numeric path: %.3e ; AD vs EXACT: %.3e\n",
              max(abs(g_ad - g_fd)), max(abs(g_ad - g_ex))))
  cat(sprintf("     sign of d/dlk21 matches the exact limit: %s\n",
              identical(sign(g_ad[3]), sign(g_ex[3]))))
}
cat("\n  n_ss = 5, same comparison:\n")
for (ext in c(FALSE, TRUE)) {
  tp <- tape_of(ext, 5L)
  g_ad <- as.numeric(tp$jacfun()(x0))
  cat(sprintf("  %-14s %12.4f %12.4f %12.4f %12.4f  |AD-exact| %.3e\n",
              paste0("AD ext=", ext), g_ad[1], g_ad[2], g_ad[3],
              g_ad[4], max(abs(g_ad - g_ex))))
}
cat("\ndone\n")
