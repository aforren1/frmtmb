# Does the order of the factor matter on frm()'s tape? gpgap's review measured it on
# the sampling route. Here: the fit's stationarity, its agreement with the closed
# form, the tape's gradient against the closed form's at one point, the outer
# standard errors, and the standard errors of eta, at true sigma 1e-3 and 3e-4,
# sorted against farthest-first, one process per arm. Designs: gpgap's seed 6601 at
# n 15 and 40, and dev/ncgp-truth.R's low-noise design, seed 20261102, at n 60 and 200.
#   Rscript dev/ncgp-perm.R <sorted|perm>
arm <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "arm", arm, "\n")
if (arm == "perm") source("dev/ncgp-perm-patch.R")
# marginal() and closed_max() of the truth script, with nothing else from it
ex <- parse("dev/ncgp-truth.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]])[1] %in% c("marginal", "closed_max", "REL")) eval(e)
}
rel_err <- function(a, b) max(abs(a / b - 1))
run <- function(tag, x, y) {
  d <- data.frame(x = x, y = y)
  w <- character()
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    frm(bf(y ~ gp(x)), data = d, family = gaussian()),
    warning = function(cnd) {
      w <<- c(w, class(cnd)[1])
      invokeRestart("muffleWarning")
    })
  el <- as.numeric(Sys.time() - t0, units = "secs")
  u <- sort(unique(x))
  J <- match(x, u)
  D2 <- outer(u, u, "-")^2
  th <- fit$estimates$theta
  pf <- c(unname(fit$estimates$beta[1]), th[1], th[2],
          unname(fit$estimates$betad[1]))
  cm <- closed_max(y, x, extra_start = pf)
  mg <- function(p) marginal(p, y, D2, J)
  # the outer order: beta, betad, theta
  to_outer <- function(p) c(p[1], p[4], p[2], p[3])
  grf <- fit$obj$gr(fit$opt$par)
  grcf_tape <- as.numeric(fit$obj$gr(to_outer(cm$par)))
  grcf_true <- numDeriv::grad(mg, cm$par, method = "Richardson")[c(1, 4, 2, 3)]
  fncf <- fit$obj$fn(to_outer(cm$par))
  # the objective's roughness along the length scale, read off second differences
  pp <- to_outer(cm$par)
  fnv <- vapply(seq_len(21L) - 11L, function(i) {
    fit$obj$fn(pp + c(0, 0, 0, i * 1e-6))
  }, 0)
  rough <- stats::sd(diff(fnv, differences = 3L)) / sqrt(20)
  H <- numDeriv::hessian(mg, cm$par)
  S <- solve(H)
  secf <- sqrt(diag(S))[c(1, 4, 2, 3)]
  sdr <- suppressWarnings(sqrt(diag(stats::vcov(fit, full = TRUE))))
  # eta's standard error at new positions: kriging variance at the closed form's
  # point plus the delta term over its parameters
  nd <- data.frame(x = seq(min(x) + 0.13, max(x) - 0.17, length.out = 7L))
  cf_eta <- function(p) {
    K <- exp(2 * p[2]) * exp(-D2 / (2 * exp(2 * p[3])))
    ks <- exp(2 * p[2]) * exp(-outer(nd$x, u, "-")^2 / (2 * exp(2 * p[3])))
    V <- K[J, J, drop = FALSE]
    diag(V) <- diag(V) * (1 + 1e-12) + exp(2 * p[4])
    R <- chol(V)
    W <- backsolve(R, t(ks[, J, drop = FALSE]), transpose = TRUE)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    list(m = p[1] + as.vector(t(W) %*% r),
         v = exp(2 * p[2]) - colSums(W^2))
  }
  ce <- cf_eta(cm$par)
  G <- numDeriv::jacobian(function(p) cf_eta(p)$m, cm$par)
  secfeta <- sqrt(pmax(ce$v, 0) + rowSums((G %*% S) * G))
  fv <- suppressWarnings(stats::fitted(fit, newdata = nd, summary = TRUE))
  cat(sprintf(paste0("NCGPO %-8s %-22s conv %d %5.2fs iter %3d evals %3d ",
                     "| gap %9.2e | max|gr| fit %8.2e, at cf tape %8.2e, true %8.2e, err %8.2e ",
                     "| rough %8.2e | fn at cf - cf %9.2e ",
                     "| outer se / cf: %s | eta se / cf: %.4f %.4f | eta est/cf sd %.2e | sigma %.4g | warn %s\n"),
              arm, tag, fit$opt$convergence, el, fit$opt$iterations,
              fit$opt$evaluations[["function"]],
              -fit$opt$objective + cm$objective, max(abs(grf)),
              max(abs(grcf_tape)), max(abs(grcf_true)),
              max(abs(grcf_tape - grcf_true)), rough, fncf - cm$objective,
              paste(sprintf("%.4f", as.numeric(sdr[seq_along(secf)]) / secf), collapse = " "),
              min(fv[, 2L] / secfeta), max(fv[, 2L] / secfeta),
              max(abs(fv[, 1L] - ce$m) / secfeta),
              exp(fit$estimates$betad[1]),
              if (length(w)) paste(unique(w), collapse = ",") else "-"))
}
for (n in c(15L, 40L)) {
  set.seed(6601)
  x <- round(stats::runif(n, 0, 6), 3)
  e0 <- stats::rnorm(n)
  for (s in c(1e-3, 3e-4)) {
    run(sprintf("seed6601 n%d s=%.0e", n, s), x, sin(2 * x) + s * e0)
  }
}
for (n in c(60L, 200L)) {
  set.seed(20261102)
  x <- round(stats::runif(n, 0, 6), 3)
  e0 <- stats::rnorm(n)
  for (s in c(1e-3, 3e-4)) {
    run(sprintf("seed2026 n%d s=%.0e", n, s), x, sin(2 * x) + s * e0)
  }
}
