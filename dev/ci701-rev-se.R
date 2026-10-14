# Review of lane ci701: standard errors of fits whose us correlation runs
# to +/-1, against an independent rank-one reference (closed-form
# Gaussian marginal likelihood, numDeriv Richardson Hessian; no RTMB).
# Prints, per case: t, lost set, frmtmb's shown SEs (outer), the
# conditional SEs from the fit-point Hessian with the correlation row
# dropped, the same at the clamp, and the reference.
#
#   Rscript dev/ci701-rev-se.R <lane|base>
args <- commandArgs(TRUE)
arm <- args[1]
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lane <- "C:/Users/adf44/source/r/wt-ci701-lib"
.libPaths(if (arm == "base") c(base, ul) else c(lane, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")

ref_fit <- function(nll, p0) {
  o <- stats::nlminb(p0, nll, control = list(eval.max = 5000,
                                              iter.max = 5000,
                                              rel.tol = 1e-14))
  o2 <- stats::optim(o$par, nll, method = "BFGS",
                     control = list(reltol = 1e-15, maxit = 5000))
  p <- if (o2$value < o$objective) o2$par else o$par
  H <- numDeriv::hessian(nll, p)
  list(par = p, nll = nll(p), se = sqrt(diag(solve(H))))
}

# rank-one mvbf (1 | q | id): stacked per id, D + z z'
v29_data <- function(sgn = 1) {
  set.seed(61)
  ng <- 12
  u <- stats::rnorm(ng, 0, 0.8)
  per <- 5
  dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                   x = stats::rnorm(ng * per))
  dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.4)
  dd$y2 <- -1 + 0.2 * dd$x + sgn * 0.6 * u[as.integer(dd$id)] +
    stats::rnorm(nrow(dd), 0, 0.3)
  dd
}
v29_ref <- function(dd, sgn) {
  gi <- split(seq_len(nrow(dd)), dd$id)
  nll <- function(p) {
    a <- p[1:2]; b <- p[3:4]; s <- exp(p[5:6]) * c(1, sgn)
    sg <- exp(p[7:8])
    tot <- 0
    for (ix in gi) {
      r <- c(dd$y1[ix] - a[1] - b[1] * dd$x[ix],
             dd$y2[ix] - a[2] - b[2] * dd$x[ix])
      n <- length(ix)
      dinv <- rep(1 / sg^2, each = n)
      z <- rep(s, each = n)
      q <- sum(z^2 * dinv)
      zr <- sum(z * dinv * r)
      ld <- sum(log(1 / dinv)) + log1p(q)
      quad <- sum(r^2 * dinv) - zr^2 / (1 + q)
      tot <- tot + 0.5 * (ld + quad + 2 * n * log(2 * pi))
    }
    tot
  }
  ref_fit(nll, c(0, 0, 0, 0, -0.3, -0.8, -1, -1))
}

rs20_data <- function(s) {
  set.seed(s)
  ng <- 20
  g <- factor(rep(seq_len(ng), each = 10))
  d <- data.frame(g = g, x = stats::rnorm(200))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(ng, 0, 0.7)[g] + stats::rnorm(200)
  d
}
rs20_ref <- function(d, sgn) {
  gi <- split(seq_len(nrow(d)), d$g)
  nll <- function(p) {
    b <- p[1:2]; s0 <- exp(p[3]); s1 <- sgn * exp(p[4]); sg <- exp(p[5])
    tot <- 0
    for (ix in gi) {
      r <- d$y[ix] - b[1] - b[2] * d$x[ix]
      z <- s0 + s1 * d$x[ix]
      n <- length(ix)
      q <- sum(z^2) / sg^2
      zr <- sum(z * r) / sg^2
      tot <- tot + 0.5 * (2 * n * log(sg) + log1p(q) +
                            sum(r^2) / sg^2 - zr^2 / (1 + q) +
                            n * log(2 * pi))
    }
    tot
  }
  ref_fit(nll, c(1, 0.5, log(0.7), -3, 0))
}

cond_se <- function(fit, j, pull) {
  h <- if (pull) {
    p <- ns$sdr_outer_point(fit$obj)
    p[j] <- sign(p[j]) * 0.9999 / sqrt(1 - 0.9999^2)
    ns$fd_hessian(p, fit$obj$gr)$H
  } else {
    ns$fd_hessian(ns$sdr_outer_point(fit$obj), fit$obj$gr)$H
  }
  sqrt(diag(solve(h[-j, -j])))
}
show <- function(tag, fit, r) {
  nm <- names(fit$opt$par)
  th <- which(nm == "theta")
  j <- th[length(th)]
  # match the reference's parameters to frmtmb's by estimate
  pf <- fit$opt$par[-j]
  ix <- vapply(pf, function(v) which.min(abs(r$par - v)), 1L)
  cat(sprintf("  match: maxdiff=%.3g dup=%d\n", max(abs(r$par[ix] - pf)),
              anyDuplicated(ix)))
  ref <- r$se[ix]
  lost <- ns$sdr_of(fit)$se_lost
  sh <- suppressWarnings(sqrt(diag(ns$sdr_of(fit)$cov.fixed)))
  cf <- cond_se(fit, j, FALSE)
  cc <- cond_se(fit, j, TRUE)
  cat(sprintf("RES %s code=%d t=%.5g lost=%s\n", tag, fit$opt$convergence,
              fit$opt$par[j],
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none"))
  cat("  names :", paste(nm[-j], collapse = " "), "\n")
  cat("  shown :", sprintf("%.6f", sh[-j]), "\n")
  cat("  condfit:", sprintf("%.6f", cf), "\n")
  cat("  condclp:", sprintf("%.6f", cc), "\n")
  cat("  ref   :", sprintf("%.6f", ref), "\n")
  cat(sprintf("  maxrel shown/ref-1=%.3g condfit=%.3g condclp=%.3g\n",
              max(abs(sh[-j] / ref - 1)), max(abs(cf / ref - 1)),
              max(abs(cc / ref - 1))))
}

f2 <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
           bf(y2 ~ x + (1 | q | id)) + gaussian())
for (sgn in c(1, -1)) {
  dd <- v29_data(sgn)
  r <- v29_ref(dd, sgn)
  cat(sprintf("REF v29 sgn=%d nll=%.10f\n", sgn, r$nll))
  for (t0 in list(NULL, 500 * sgn)) {
    st <- if (!is.null(t0)) list(theta = c(0, 0, t0))
    fit <- suppressMessages(suppressWarnings(frm(f2, data = dd,
                                                 start = st)))
    cat(sprintf("  frmtmb logLik=%.10f\n", as.numeric(logLik(fit))))
    show(sprintf("v29 sgn=%d t0=%s", sgn, if (is.null(t0)) "def" else t0),
         fit, r)
  }
}
for (s in c(2L, 7L, 8L, 11L, 12L, 15L, 17L, 19L)) {
  d <- rs20_data(s)
  fit <- suppressMessages(suppressWarnings(
    frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = d)))
  sgn <- sign(fit$estimates$theta[3])
  r <- rs20_ref(d, sgn)
  cat(sprintf("REF rs20 seed=%d sgn=%d nll=%.10f frm=%.10f\n", s, sgn,
              r$nll, -as.numeric(logLik(fit))))
  show(sprintf("rs20 seed=%d", s), fit, r)
}
