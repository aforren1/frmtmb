# Two questions from dev/ncgp-perm.R and dev/ncgp-perm2.R. (1) perm.R read eta's
# standard error 12.6 times the closed form's on sorted positions at one design and
# perm2.R reads 1.02 at most on thirty others: is 12.6 the fit or the reference? The
# reference read at frm()'s point and at the closed form's maximum, with the worst
# position. (2) conv 1 at true sigma 3e-4 on half of perm2's rows: how far short, and
# what warning does the user see? Same designs on every build.
#   Rscript dev/ncgp-perm3.R <sorted|perm|base>
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm != "base") libs <- c("C:/Users/adf44/source/r/wt-ncgp-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "arm", arm, "\n")
if (arm == "perm") source("dev/ncgp-perm-patch.R")
ex <- parse("dev/ncgp-truth.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]])[1] %in% c("marginal", "closed_max", "REL")) eval(e)
}
one <- function(tag, x, y, m_grid) {
  d <- data.frame(x = x, y = y)
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ gp(x)), data = d, family = gaussian()),
    warning = function(cnd) {
      w <<- c(w, substr(conditionMessage(cnd), 1, 60))
      invokeRestart("muffleWarning")
    })
  u <- sort(unique(x))
  J <- match(x, u)
  D2 <- outer(u, u, "-")^2
  th <- fit$estimates$theta
  pf <- c(unname(fit$estimates$beta[1]), th[1], th[2], unname(fit$estimates$betad[1]))
  cm <- closed_max(y, x, extra_start = pf)
  # the base's objective has another nugget, so the gap is read on one marginal
  gap <- -marginal(pf, y, D2, J) + cm$objective
  nd <- data.frame(x = seq(min(x) + 0.13, max(x) - 0.17, length.out = m_grid))
  ref_at <- function(p0) {
    S <- solve(numDeriv::hessian(function(p) marginal(p, y, D2, J), p0))
    cf <- function(p) {
      K <- exp(2 * p[2]) * exp(-D2 / (2 * exp(2 * p[3])))
      ks <- exp(2 * p[2]) * exp(-outer(nd$x, u, "-")^2 / (2 * exp(2 * p[3])))
      V <- K[J, J, drop = FALSE]
      diag(V) <- diag(V) * (1 + 1e-12) + exp(2 * p[4])
      R <- chol(V)
      W <- backsolve(R, t(ks[, J, drop = FALSE]), transpose = TRUE)
      list(m = p[1] + as.vector(t(W) %*% backsolve(R, y - p[1], transpose = TRUE)),
           v = exp(2 * p[2]) - colSums(W^2))
    }
    G <- numDeriv::jacobian(function(p) cf(p)$m, p0)
    sqrt(pmax(cf(p0)$v, 0) + rowSums((G %*% S) * G))
  }
  fv <- suppressWarnings(stats::fitted(fit, newdata = nd))
  rf <- ref_at(pf)
  rc <- ref_at(cm$par)
  r1 <- fv[, 2L] / rf
  r2 <- fv[, 2L] / rc
  k <- which.max(abs(r1 - 1))
  cat(sprintf(paste0("NCGPR3 %-6s %-20s conv %d ll gap %9.2e max|gr| %8.2e sigma %.4g cf %.4g | ",
                     "se/ref(frm) %.4f %.4f at x %.3f (se %.3e ref %.3e) | se/ref(cf) %.4f %.4f | ",
                     "warn %s\n"),
              arm, tag, fit$opt$convergence, gap, max(abs(fit$obj$gr(fit$opt$par))),
              exp(fit$estimates$betad[1]), exp(cm$par[4]),
              min(r1), max(r1), nd$x[k], fv[k, 2L], rf[k], min(r2), max(r2),
              if (length(w)) paste(unique(w), collapse = "; ") else "-"))
}
# perm.R's rows, its 7-point grid
for (n in c(40L)) {
  set.seed(6601)
  x <- round(stats::runif(n, 0, 6), 3)
  e0 <- stats::rnorm(n)
  one(sprintf("seed6601 n%d 3e-4", n), x, sin(2 * x) + 3e-4 * e0, 7L)
}
for (n in c(200L)) {
  set.seed(20261102)
  x <- round(stats::runif(n, 0, 6), 3)
  e0 <- stats::rnorm(n)
  for (s in c(1e-3, 3e-4)) {
    one(sprintf("seed2026 n%d %.0e", n, s), x, sin(2 * x) + s * e0, 7L)
  }
}
# perm2's designs, order of drawing as there, the conv question
for (n in c(40L, 100L, 200L)) for (seed in 1:5) {
  set.seed(seed)
  x <- round(stats::runif(n, 0, 6), 3)
  y <- sin(2 * x) + 3e-4 * stats::rnorm(n)
  one(sprintf("n%d seed%d 3e-4", n, seed), x, y, 25L)
}
