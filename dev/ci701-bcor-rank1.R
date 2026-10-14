# Lane ci701: are the two sds of the (1 | q | id) block identified when
# the correlation is at +1? The model at the boundary is a rank-one
# random effect, y1 = a1 + b1 x + s1 u + e1, y2 = a2 + b2 x + s2 u + e2,
# u ~ N(0, 1) per id, which RTMB fits with no near-singular covariance.
# Its maximum must equal frmtmb's log-likelihood at the boundary, and its
# Hessian gives the sds' standard errors at the boundary. Then frmtmb's
# own objective is evaluated with the correlation parameter pulled back
# to t = 10, 30, 70.7 (r = 0.9999), 300, 1000, 3000: the FD Hessian of the
# two log sds, the gradient's spread over ulp-moved points.
#
#   Rscript dev/ci701-bcor-rank1.R <lib|base>
args <- commandArgs(TRUE)
lib <- args[1]
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
library(RTMB)
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

gi <- as.integer(dd$id)
pars <- list(a = c(0, 0), b = c(0, 0), ls = c(0, 0), le = c(0, 0),
             u = numeric(ng))
nll <- function(p) {
  getAll(p)
  m1 <- a[1] + b[1] * dd$x + exp(ls[1]) * u[gi]
  m2 <- a[2] + b[2] * dd$x + exp(ls[2]) * u[gi]
  -sum(dnorm(dd$y1, m1, exp(le[1]), log = TRUE)) -
    sum(dnorm(dd$y2, m2, exp(le[2]), log = TRUE)) -
    sum(dnorm(u, 0, 1, log = TRUE))
}
obj <- MakeADFun(nll, pars, random = "u", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
cat(sprintf("RANK1 code=%d msg=%s ll=%.10f\n", opt$convergence, opt$message,
            -opt$objective))
s <- summary(sdr, "fixed")
print(s)
cat(sprintf("RANK1 log sd1 %.6f se %.6f; log sd2 %.6f se %.6f\n",
            s[5, 1], s[5, 2], s[6, 1], s[6, 2]))

f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
fit <- suppressWarnings(frm(f, data = dd))
cat(sprintf("FRM code=%d ll=%.10f theta=%s\n", fit$opt$convergence,
            as.numeric(logLik(fit)),
            paste(sprintf("%.6f", fit$estimates$theta), collapse = ",")))
nm <- names(fit$opt$par)
th <- which(nm == "theta")
p0 <- fit$opt$par
for (t3 in c(10, 30, 70.7, 300, 1000, 3000, 10000)) {
  p <- p0
  p[th[3]] <- t3
  # re-maximize the rest with the correlation parameter held
  fr <- setdiff(seq_along(p), th[3])
  o <- nlminb(p[fr], function(q) { p[fr] <- q; fit$obj$fn(p) },
              function(q) { p[fr] <- q; fit$obj$gr(p)[fr] })
  p[fr] <- o$par
  h <- frmtmb:::fd_hessian(p, fit$obj$gr)
  H <- h$H
  V <- tryCatch(solve(H[fr, fr]), error = function(e) NULL)
  se <- if (is.null(V)) c(NA, NA) else sqrt(diag(V))[match(th[1:2], fr)]
  gs <- vapply(1:6, function(j) {
    q <- p * (1 + sample(0:4, length(p), TRUE) * 2^-52)
    max(abs(fit$obj$gr(q)[th[1:2]]))
  }, 0)
  cat(sprintf(
    "HOLD t3=%g ll=%.10f code=%d H11=%.4g H22=%.4g H12=%.4g H33=%.3g se(th1,th2)=%.5f,%.5f max|g_th12| over 6 ulp moves=%.3g\n",
    t3, -o$objective, o$convergence, H[th[1], th[1]], H[th[2], th[2]],
    H[th[1], th[2]], H[th[3], th[3]], se[1], se[2], max(gs)))
}
