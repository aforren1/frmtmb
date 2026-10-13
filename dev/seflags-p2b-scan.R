# Lane seflags, punch round 2b: the "illcond" directions put back into
# the covariance. For raw polynomials of degree 8 on [1, 2], with and
# without (1 | g), seeds 1 to 12, and the B5 design (z and a factor w
# added): the coefficient SEs against the exact ones (QR of the whitened
# design at the fit's sigma and sd). On a build before the change the
# stored record is used, V + d diag(1 / lambda) d'; after it, vcov().
# Also the measurement that would decide trust: each direction's
# curvature by se_obj_curv() at the step it took and at twice and half
# that step.
#   Rscript dev/seflags-p2b-scan.R <lane|op3|m9|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
exact_se <- function(f, X, d, re) {
  on <- ns$outer_par_names(f)
  sg <- exp(f$opt$par[grepl("^sigma", on)])
  tg <- if (re) exp(f$opt$par[on == "theta_1"]) else 0
  Vv <- sg^2 * diag(nrow(X))
  if (re) Vv <- Vv + tg^2 * tcrossprod(stats::model.matrix(~ 0 + g, d))
  q <- qr(backsolve(chol(Vv), X, transpose = TRUE), tol = 1e-30)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  ex <- numeric(ncol(X))
  ex[q$pivot] <- sqrt(rowSums(Ri^2))
  ex
}
for (des in c("poly", "b5")) for (re in c(TRUE, FALSE)) for (s in 1:12) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  if (des == "b5") {
    d$z <- stats::rnorm(200)
    d$w <- factor(sample(c("a", "b"), 200, TRUE))
    d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
      stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  } else {
    d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
      stats::rnorm(200, 0, 0.2)
  }
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:8, collapse = " + ")
  if (des == "b5") rhs <- paste(rhs, "+ z + w")
  r <- cap(frm(stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)")),
               data = d))
  f <- r$v
  if (inherits(f, "rev_err")) {
    cat(sprintf("%s re %s seed %2d | error %s\n", des, re, s, f))
    next
  }
  sdr <- suppressWarnings(ns$sdr_of(f))
  il <- sdr$se_null_flat$ill
  on <- ns$outer_par_names(f)
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  cn <- colnames(X)
  ex <- exact_se(f, X, d, re)
  how <- "vcov"
  if (!is.null(il) && exists("se_ill_share", ns)) {
    Vf <- il$V + il$d %*% (t(il$d) / il$lam)
    se <- sqrt(diag(Vf))[match(cn, on[il$rows])]
    how <- "back"
  } else {
    V <- suppressWarnings(as.matrix(vcov(f, full = TRUE)))
    se <- sqrt(diag(V))[match(cn, on)]
  }
  lost <- sdr$se_lost
  cat(sprintf(paste0("%s re %s seed %2d | code %d | warn %d | lost %s | ",
                     "ill %s | %s SE/exact %s | %s\n"),
              des, re, s, f$opt$convergence, length(r$w),
              if (length(lost)) paste(unique(lost), collapse = ",") else "-",
              if (is.null(il)) "-" else paste(signif(il$lam, 3),
                                               collapse = " "),
              how, paste(range(signif(se / ex, 4)), collapse = " to "),
              substr(paste(r$w, collapse = " + "), 1, 45)))
  if (!is.null(il)) {
    p <- f$opt$par
    for (k in seq_along(il$lam)) {
      dn <- numeric(length(p))
      dn[il$rows] <- il$d[, k]
      c1 <- ns$se_obj_curv(f, p, dn)
      st <- attr(c1, "step")
      cc <- function(m) {
        h <- st * m
        fn <- function(x) f$obj$fn(x)
        (fn(p + h * dn) + fn(p - h * dn) - 2 * fn(p)) / h^2
      }
      # what a Newton step along the direction would gain
      gd <- sum(as.numeric(f$obj$gr(p)) * dn)
      cat(sprintf(paste0("    dir %d lam %.4g obj %.4g step %.3g | at 2x %.4g,",
                         " 0.5x %.4g, 0.25x %.4g | grad.d %.3g gain %.4g\n"),
                  k, il$lam[k], c1, st, cc(2), cc(0.5), cc(0.25), gd,
                  gd^2 / (2 * il$lam[k])))
    }
  }
}
