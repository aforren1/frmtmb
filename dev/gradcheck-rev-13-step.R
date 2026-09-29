## Reviewer, claim 1(b) and the section 5 dispute: on the collinear
## ridge at eps = 1e-6 and 1e-7 a HARD re-optimization recovers only
## 1.7e-08 log-likelihood units, while the headroom reads 0.529 and
## 0.737. Either the headroom is 3e7 times too large, or the optimizer
## cannot find a point the Newton step reaches. Settled by taking the
## step and evaluating the objective, which needs no optimizer.
## usage: Rscript gradcheck-rev-13-step.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
for (eps in c(1e-5, 1e-6, 1e-7)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  dd <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, eps))
  dd$y <- rnorm(n, 1 + 2 * dd$x1, 1)
  f <- suppressWarnings(frm(bf(y ~ x1 + x2), family = gaussian(), data = dd,
                            control = frmtmb_control(restarts = 0,
                                                     autoscale = FALSE)))
  p <- f$opt$par
  g <- drop(f$obj$gr(p))
  H <- f$obj$he(p); H <- (H + t(H)) / 2
  ch <- chol(H)
  step <- -backsolve(ch, backsolve(ch, g, transpose = TRUE))
  pred <- 0.5 * sum(backsolve(ch, g, transpose = TRUE)^2)
  f0 <- as.numeric(f$obj$fn(p))
  best <- -Inf; bt <- NA_real_
  for (t in c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0.01)) {
    drop_t <- f0 - as.numeric(f$obj$fn(p + t * step))
    if (is.finite(drop_t) && drop_t > best) { best <- drop_t; bt <- t }
  }
  cat(sprintf(paste0("eps=%-8s rcond(H)=%-11s |step|=%-11s predicted=%-11s",
                     " best drop=%-11s at t=%-5s ratio=%s\n"),
              format(eps), format(1 / kappa(H, exact = TRUE), digits = 4),
              format(sqrt(sum(step^2)), digits = 4),
              format(pred, digits = 6), format(best, digits = 6),
              format(bt), format(best / pred, digits = 4)))
  ## and is the point the step reaches a legitimate one? print the
  ## coefficients there and the objective at both points at full
  ## precision, so a catastrophic cancellation would show
  cat("   par at optimum  :", paste(format(p, digits = 10),
                                    collapse = " "), "\n")
  cat("   par after step  :", paste(format(p + bt * step, digits = 10),
                                    collapse = " "), "\n")
  cat("   objective       :", format(f0, digits = 15), " ->",
      format(as.numeric(f$obj$fn(p + bt * step)), digits = 15), "\n")
  ## the same two points scored by an independent route: lm()'s residual
  ## sum of squares under the exact gaussian log likelihood
  X <- cbind(1, dd$x1, dd$x2)
  nll <- function(q) {
    mu <- X %*% q[1:3]
    s <- exp(q[[4]])
    -sum(stats::dnorm(dd$y, mu, s, log = TRUE))
  }
  cat("   independent nll :", format(nll(p), digits = 15), " ->",
      format(nll(p + bt * step), digits = 15),
      " drop", format(nll(p) - nll(p + bt * step), digits = 6), "\n\n")
}
cat("DONE step\n")
