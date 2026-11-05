# Cost of one exact gp() fit, centered (rellib-r16) against non-centered (lane).
# y = sin(x) + 0.3 e, seed 20261103. n rows, at most 500 positions: n = 1000 has
# 500 positions, two rows each. Prints per-evaluation minima over blocks of at
# least 1.2 s, tape nodes, outer iterations, inner solves, the whole frm() and
# summary(), with a fixed arithmetic control for the machine's load. Run under
# dev/ncgp-guard.ps1, which reads the peak working set.
#   Rscript dev/ncgp-perf.R <label> <lib or "base"> <n> [family] [rounds]
args <- commandArgs(TRUE)
label <- args[1]
libs <- c("C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (args[2] != "base") libs <- c(args[2], libs)
.libPaths(libs)
n <- as.integer(args[3])
family <- if (length(args) > 3) args[4] else "gaussian"
rounds <- if (length(args) > 4) as.integer(args[5]) else 3L
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261103)
npos <- min(n, 500L)
pos <- sort(runif(npos, 0, 10))
x <- rep(pos, length.out = n)
d <- switch(
  family,
  gaussian = data.frame(x = x, y = sin(x) + 0.3 * rnorm(n)),
  poisson = data.frame(x = x, y = rpois(n, exp(1 + sin(x)))),
  bernoulli = data.frame(x = x, y = rbinom(n, 1, plogis(2 * sin(x)))))
fam <- switch(family, gaussian = gaussian(), poisson = poisson(),
              bernoulli = bernoulli())
control_block <- function() {
  m <- matrix(rnorm(1e6), 1000)
  t0 <- proc.time()[[3]]; for (i in 1:5) crossprod(m); (proc.time()[[3]] - t0) / 5
}
minblock <- function(f) {
  out <- Inf
  for (r in seq_len(rounds)) {
    t0 <- Sys.time(); k <- 0L
    repeat {
      f(); k <- k + 1L
      el <- as.numeric(Sys.time() - t0, units = "secs")
      if (el > 1.2) break
    }
    out <- min(out, el / k)
  }
  out
}
c0 <- control_block()
t0 <- Sys.time()
fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = fam))
tfit <- as.numeric(Sys.time() - t0, units = "secs")
obj <- fit$obj
nodes <- function(nm) {
  tp <- tryCatch(RTMB::GetTape(obj, nm), error = function(e) NULL)
  if (is.null(tp)) NA_real_ else nrow(tp$data.frame())
}
nn <- c(F = nodes("ADFun"), G = nodes("ADGrad"), H = nodes("ADHess"))
p0 <- fit$opt$par
i <- 0L
wob <- function() {
  i <<- i + 1L
  p0 + 1e-3 * sin(i) * seq_along(p0) / length(p0)
}
# inner solves per outer evaluation: count of spHess calls
env <- obj$env
sph <- env$spHess
cnt <- 0L
env$spHess <- function(...) { cnt <<- cnt + 1L; sph(...) }
obj$fn(wob()); c1 <- cnt
env$spHess <- sph
tfn <- minblock(function() obj$fn(wob()))
tgr <- minblock(function() obj$gr(wob()))
t1 <- Sys.time()
s <- sdr <- tryCatch(frmtmb:::autoscale_sdreport(fit, jp = TRUE), error = function(e) NULL)
tsd <- as.numeric(Sys.time() - t1, units = "secs")
t1 <- Sys.time()
fv <- suppressWarnings(fitted(fit))
tfitted <- as.numeric(Sys.time() - t1, units = "secs")
t1 <- Sys.time()
nd <- data.frame(x = seq(0.1, 9.9, length.out = 20))
pv <- suppressWarnings(frm_linpred(fit, newdata = nd, se.fit = TRUE))
tpred <- as.numeric(Sys.time() - t1, units = "secs")
c1b <- control_block()
cat(sprintf(paste0("NCGPF %s n=%d pos=%d %s fit=%.2fs iter=%d feval=%d geval=%d ",
                   "ll=%.4f sdgp=%.4f nodes(F,G,H)=%s fn=%.4fs gr=%.4fs ",
                   "hessPerEval=%d sdreport=%.2fs fitted=%.2fs predict=%.2fs ",
                   "control=%.3f,%.3f\n"),
            label, n, npos, family, tfit, fit$opt$iterations,
            fit$opt$evaluations[[1]], fit$opt$evaluations[[2]],
            -fit$opt$objective,
            exp(fit$estimates$theta[1]), paste(nn, collapse = ","),
            tfn, tgr, c1, tsd, tfitted, tpred, c0, c1b))
