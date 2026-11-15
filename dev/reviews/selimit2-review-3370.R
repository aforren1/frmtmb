# Review of lane selimit2, item 3370 (joint_cov_white_inner()).
#   Rscript selimit2-review-3370.R <base|lane> exact
#     fitted(re_formula = NA) standard errors of a REML fit on a raw
#     polynomial (degrees 5, 6, 7) with (1 | g) and no group effect, against
#     a (X' W X)^-1 a from the weighted QR of an orthogonal basis, for gaussian()
#     and for gaussian() with sigma ~ f (weights 1 / sigma_i^2)
#   Rscript selimit2-review-3370.R <base|lane> cost <groups>
#     time and peak memory of the first fitted(re_formula = NA) on a large
#     fit of the same kind
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
mode <- a[2]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    cat("   W:", substr(conditionMessage(w), 1, 120), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
}
mk <- function(n, ng, deg, seed, sig2 = FALSE) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(seq_len(ng), length.out = n)),
                  f = factor(rep(c("a", "b"), length.out = n)))
  s <- if (sig2) ifelse(d$f == "a", 0.2, 0.5) else 0.3
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, s)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  d
}
if (mode == "exact") {

  for (sig2 in c(FALSE, TRUE)) {
    for (deg in 5:7) {
      for (seed in 1:3) {
        d <- mk(200, 20, deg, seed, sig2)
        rhs <- paste(paste0("x", 1:deg), collapse = " + ")
        fo <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
        bfo <- if (sig2) bf(fo, sigma ~ f) else bf(fo)
        fit <- quiet(frm(bfo, data = d, REML = TRUE))
        lost <- frmtmb:::sdr_of(fit)$se_lost
        th <- fit$opt$par
        X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
        nd <- data.frame(x = seq(1, 2, by = 0.25), g = d$g[1], f = "a")
        for (k in 1:deg) nd[[paste0("x", k)]] <- nd$x^k
        A <- stats::model.matrix(stats::as.formula(paste("~", rhs)), nd)
        # sigma per row at the estimate (sigma ~ f or a constant)
        sg <- if (sig2) {
          bd <- fit$estimates$betad
          exp(bd[1] + bd[2] * (d$f == "b"))
        } else rep(exp(th[["betad"]]), nrow(d))
        # the group sd at its estimate, as a variance ratio per row's sigma:
        # with it lost at ~0 the GLS weights are 1 / sigma_i^2
        # the same column space in an orthogonal polynomial basis, whose
        # weighted QR is well conditioned: a' (X' W X)^-1 a is invariant
        pb <- stats::poly(d$x, deg)
        Xo <- cbind(1, pb) / sg
        Ao <- cbind(1, stats::predict(pb, nd$x))
        Ri <- backsolve(qr.R(qr(Xo)), diag(ncol(Xo)))
        vref <- rowSums((Ao %*% Ri)^2)
        if (sig2) {
          # what the sigma parameters' uncertainty adds through the GLS
          # weights: (A dbeta/ds) Vs (A dbeta/ds)'
          pred <- function(b) {
            s <- exp(b[1] + b[2] * (d$f == "b"))
            q <- qr(cbind(1, pb) / s)
            drop(Ao %*% qr.coef(q, d$y / s))
          }
          bd <- unname(fit$estimates$betad)
          h <- 1e-5
          Jp <- vapply(1:2, function(k) {
            e <- replace(numeric(2), k, h)
            (pred(bd + e) - pred(bd - e)) / (2 * h)
          }, numeric(nrow(Ao)))
          Vf <- quiet(stats::vcov(fit, full = TRUE))
          ib <- grep("^sigma_", rownames(Vf))
          Vs <- Vf[ib, ib]
          vref <- vref + rowSums((Jp %*% Vs) * Jp)
        }
        fp <-quiet(stats::fitted(fit, newdata = nd, re_formula = NA))
        r <- unname(fp[, "Est.Error"]) / sqrt(vref) - 1
        cat(sprintf(paste0("%s deg %d seed %d | sd est %.3g | lost %s | ",
                           "whiten %s | max |SE/ref - 1| %.3g\n"),
                    if (sig2) "sigma ~ f" else "gaussian", deg, seed,
                    exp(fit$estimates$theta[1]),
                    paste(names(lost), lost, collapse = ","),
                    !is.null(fit$beta_whiten), max(abs(r))))
      }
    }
  }
} else {
  ng <- as.integer(a[3])
  d <- mk(2L * ng, ng, 7, 1)
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g)
  t0 <- proc.time()[["elapsed"]]
  fit <- quiet(frm(fo, data = d, REML = TRUE))
  cat("fit s", round(proc.time()[["elapsed"]] - t0, 2), "| whiten",
      !is.null(fit$beta_whiten), "| lost",
      paste(names(frmtmb:::sdr_of(fit)$se_lost),
            frmtmb:::sdr_of(fit)$se_lost), "\n")
  nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
  for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  invisible(gc(reset = TRUE))
  t0 <- proc.time()[["elapsed"]]
  fp <- quiet(stats::fitted(fit, newdata = nd, re_formula = NA))
  el <- proc.time()[["elapsed"]] - t0
  mx <- gc()[, 6]
  cat(sprintf("groups %d | joint size %d | fitted() s %.2f | max used MB %.0f\n",
              ng, nrow(frmtmb:::get_joint_cov(fit)$V), el, sum(mx)))
  print(fp[, c("Estimate", "Est.Error")])
}
