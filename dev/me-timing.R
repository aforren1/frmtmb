# Does the me() change add work to the objective of a model WITHOUT
# me() terms? (lane me, 2026-09-25)
#
# The objective gains two things for every model: a NULL check on
# frame[["me"]] per evaluation and an empty loop over lp[["me"]] in
# lp_eta_fixed(). Neither puts a node on the tape, so the load-independent
# count is the tape itself: the number of operations RTMB records, and
# the value and gradient at a fixed parameter vector, which must be
# bitwise identical between the base build and the lane build. The
# clock is reported beside it, as a minimum over rounds, with the same
# arithmetic control run in both processes.
#
#   Rscript dev/me-timing.R base   > dev/me-timing-base.txt
#   Rscript dev/me-timing.R lane   > dev/me-timing-lane.txt
lib <- if (identical(commandArgs(TRUE)[1], "base")) character(0) else
  "/opt/rlib/lane-me"
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")

set.seed(2026)
n <- 3000
ng <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(seq_len(ng), n, TRUE)))
u <- matrix(rnorm(2 * ng, 0, c(0.8, 0.3)), ng, 2, byrow = TRUE)
d$y <- 1 + 0.5 * d$x - 0.3 * d$z + u[d$g, 1] + u[d$g, 2] * d$x +
  rnorm(n, 0, exp(0.2 * d$z))
form <- bf(y ~ x + z + (1 + x | g), sigma ~ z) + gaussian()

fit <- frm(form, data = d)
tape <- RTMB::GetTape(fit$obj)
ops <- tryCatch(length(tape$get_graph()$from),
                error = function(e) NA_integer_)
cat("tape print:\n")
print(tape)
p <- fit$obj$env$last.par.best
v <- fit$obj$env$f(p)
g <- fit$obj$env$f(p, order = 1)
cat(sprintf("logLik %.17g\n", as.numeric(logLik(fit))))
cat(sprintf("f(last.par.best) %.17g\n", v))
cat(sprintf("sum|grad| %.17g\n", sum(abs(g))))
cat("ops", ops, "\n")

# clock: whole fits, minimum over rounds, with an arithmetic control
ctl <- function() {
  s <- 0
  for (i in 1:2e6) s <- s + sqrt(i)
  s
}
tf <- tc <- numeric(0)
for (r in 1:5) {
  t0 <- proc.time()[["elapsed"]]
  invisible(frm(form, data = d))
  tf <- c(tf, proc.time()[["elapsed"]] - t0)
  t0 <- proc.time()[["elapsed"]]
  invisible(ctl())
  tc <- c(tc, proc.time()[["elapsed"]] - t0)
}
cat(sprintf("fit seconds min %.3f (rounds %s)\n", min(tf),
            paste(sprintf("%.2f", tf), collapse = " ")))
cat(sprintf("control seconds min %.3f (rounds %s)\n", min(tc),
            paste(sprintf("%.2f", tc), collapse = " ")))
