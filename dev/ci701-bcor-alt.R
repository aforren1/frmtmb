# Lane ci701: two candidate repairs for the boundary-correlation fit of
# test-review-v29.R (set.seed(61)), measured on data moved by a few ulps
# (seed 1000 + k, as dev/ci701-bcor-probe.R):
#   a. the correlation parameter at its limit removed, the rest of the
#      finite-difference Hessian taken at the fit's point;
#   b. the same, with the Hessian taken where the correlation parameter
#      is pulled back to the clamp, t = r / sqrt(1 - r^2), r = 0.9999.
# Reference: the rank-one fit of dev/ci701-bcor-rank1.R, se(log sd1)
# 0.217708, se(log sd2) 0.222476.
#
#   Rscript dev/ci701-bcor-alt.R <lib|base> <k0> <k1>
args <- commandArgs(TRUE)
lib <- args[1]; k0 <- as.integer(args[2]); k1 <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")

set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
d0 <- dd
mdata <- function(k) {
  dd <- d0
  if (k > 0) {
    set.seed(1000 + k)
    for (v in c("x", "y1", "y2")) {
      dd[[v]] <- dd[[v]] * (1 + sample(0:4, nrow(dd), TRUE) * 2^-52)
    }
  }
  dd
}
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
r <- 0.9999
tc <- r / sqrt(1 - r^2)
ses <- function(H, th) {
  fr <- setdiff(seq_len(nrow(H)), th[3])
  V <- tryCatch(solve(H[fr, fr]), error = function(e) NULL)
  if (is.null(V)) return(c(NA, NA))
  d <- diag(V)[match(th[1:2], fr)]
  ifelse(d > 0, sqrt(d), NaN)
}
for (k in k0:k1) {
  fit <- suppressWarnings(frm(f, data = mdata(k)))
  p <- fit$opt$par
  th <- which(names(p) == "theta")
  g <- drop(fit$obj$gr(p))
  H0 <- frmtmb:::fd_hessian(p, fit$obj$gr)$H
  q <- p
  q[th[3]] <- sign(p[th[3]]) * min(abs(p[th[3]]), tc)
  H1 <- frmtmb:::fd_hessian(q, fit$obj$gr)$H
  cat(sprintf("RES k=%d code=%d gmax=%.3g t=%.1f a=%s b=%s full=%s\n", k,
              fit$opt$convergence, max(abs(g)), p[th[3]],
              paste(sprintf("%.5f", ses(H0, th)), collapse = ","),
              paste(sprintf("%.5f", ses(H1, th)), collapse = ","),
              paste(sprintf("%.5f", (function(v) ifelse(v > 0, sqrt(v), NaN))(diag(solve(H0)))[th[1:2]]),
                    collapse = ",")))
}
