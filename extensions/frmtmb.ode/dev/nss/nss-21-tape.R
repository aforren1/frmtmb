# NIT 6: the tape, which the lane claimed nothing about and should have.
#
# The node-count method is the reviewer's, from
# `dev/rev-nss/rev-nss-04-cost-grad.R` 3b: the printed operation stack
# is one line per node, so its length is the node count, and a count is
# load-independent where a clock is not.
#
# Run with NSS_ARM=ref for the base commit.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-21-tape.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
nss_report_env()
ARM <- Sys.getenv("NSS_ARM", "lane")

two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}
TT <- c(0, 6, 12, 23.9)
EV <- data.frame(time = 0, state = 1L, value = 100, ii = 24, ss = TRUE)
f_ode <- function(th, ext, n) {
  args <- list(two_oral, init = list(0, 0, 0), times = TT,
               parms = list(exp(th[1L]), exp(th[2L]), exp(th[3L]),
                            exp(th[4L])),
               events = EV, output = 2L, n_ss = n, ss_tol = Inf,
               atol = 1e-8, rtol = 1e-8)
  if (ARM != "ref") args$ss_extrapolate <- ext
  sum(do.call(frm_ode, args))
}
x0 <- log(c(0.15, 0.3, 0.02, 1.0))
n_nodes <- function(tp) length(capture.output(tp$print(depth = 1))) - 1L
tape_of <- function(ext, n) MakeTape(function(th) {
  "c" <- RTMB::ADoverload("c")
  f_ode(th, ext, n)
}, x0)

cat("\n== nodes in the outer tape ==\n")
cat(sprintf("%6s %14s %14s %10s\n", "n_ss", "truncated",
            "extrapolated", "delta"))
for (n in c(4L, 5L, 20L, 40L)) {
  a <- n_nodes(tape_of(FALSE, n))
  b <- n_nodes(tape_of(TRUE, n))
  cat(sprintf("%6d %14d %14d %+10d (%+.1f%%)\n", n, a, b, b - a,
              100 * (b - a) / a))
}

if (ARM != "ref") {
  cat("\n== one ode_ss_extrapolate() on three states ==\n")
  solo <- n_nodes(MakeTape(function(x) sum(x), rep(1, 9)))
  withc <- n_nodes(MakeTape(function(x) {
    "c" <- RTMB::ADoverload("c")
    z <- frmtmb.ode:::ode_ss_extrapolate(x[1:3], x[4:6], x[7:9],
                                         1e-8, 1e-8)
    sum(z[["y"]]) + sum(z[["r"]])
  }, rep(1, 9)))
  cat("  ", withc - solo, "nodes\n")
}
