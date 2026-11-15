# Review of lane arspat2: the fcor() margin route failed (NaN gradient) at
# n = 1000, m = 200 where 0.82.0 fit. Smaller designs: k x k rook lattice
# (row-standardized, M = solve(crossprod(I - 0.4 W))), m missing, (1 | h),
# student() with gamma(2, 0.1) on nu, data seed 11, missing seed m (as the
# cost script). Base when REVLIB="".
source("dev/reviews/arspat2-review-common.R")
cells <- list(c(20, 100), c(20, 200), c(20, 300))
for (cl in cells) {
  k <- cl[1]; m <- cl[2]
  W <- rook(k); W <- W / rowSums(W)
  n <- nrow(W)
  A0 <- diag(n) - 0.4 * W
  M <- solve(crossprod(A0))
  set.seed(11)
  d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
  d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(8)[d$h] +
                           rnorm(n) * sqrt(4 / rchisq(1, 4))))
  set.seed(m)
  d$y[sort(sample(n, m))] <- NA
  pr <- set_prior("gamma(2, 0.1)", class = "nu")
  for (ty in c("fcor", "sar")) {
    fo <- if (ty == "fcor") y | mi() ~ x + (1 | h) + fcor(M) else
      y | mi() ~ x + (1 | h) + sar(W, type = "error")
    d2 <- if (ty == "fcor") list(M = M) else list(W = W)
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(quiet(frm(fo, data = d, data2 = d2, family = student(),
                              prior = pr)), error = function(e) e)
    tt <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      cat(sprintf("n %d m %d %s: ERROR after %.1f s: %s\n", n, m, ty, tt,
                  substr(conditionMessage(fit), 1, 90)))
      # where: the objective at the start, verbose dry run
      ob <- tryCatch(quiet(frm(fo, data = d, data2 = d2, family = student(),
                               prior = pr, dry_run = "objective")),
                     error = function(e) e)
      if (!inherits(ob, "error")) {
        p0 <- ob$obj$par
        f0 <- ob$obj$fn(p0)
        g0 <- ob$obj$gr(p0)
        cat(sprintf("   at start: fn %s, gr finite %s, inner max |miss| %.3g\n",
                    format(f0), all(is.finite(g0)),
                    max(abs(ob$obj$env$last.par[ob$obj$env$random]))))
        tr <- tryCatch(quiet(frm(fo, data = d, data2 = d2, family = student(),
                                 prior = pr, control = frmtmb_control(
                                   optimizer = "optim"))),
                       error = function(e) e)
        cat("   optim:", if (inherits(tr, "error")) substr(conditionMessage(tr), 1, 80)
            else sprintf("ok logLik %.4f", as.numeric(logLik(tr))), "\n")
      }
    } else {
      cat(sprintf("n %d m %d %s: ok %.1f s logLik %.4f conv %d\n", n, m, ty, tt,
                  as.numeric(logLik(fit)), fit$opt$convergence))
    }
  }
}
