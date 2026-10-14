# Lane illse: covariance accuracy on bernoulli, poisson and cumulative
# raw polynomials (degrees 6 to 8, with and without (1 | g)).
# References:
#   exact  without random effects, the Hessian of the log-likelihood at
#          the fit's own point in 256 bits (bernoulli and poisson: X'WX,
#          canonical links; cumulative: second differences of its
#          log-likelihood in the thresholds frmtmb estimates), inverted
#          there;
#   poly   the same likelihood fitted with poly(x, deg), well conditioned,
#          its covariance mapped to the raw coefficients through the
#          exact change of basis X_raw = X_poly B (256 bits); with
#          (1 | g) this is the only reference. Its own point differs
#          from the raw fit's by the two optimizers' tolerances.
# Per fit: code, warnings, lost set, Newton steps, the log-likelihood
# gap to the poly() fit, entry errors against each reference (cor, rel,
# diag as in dev/illse-check.R, over the coefficient block), and the
# Est.Error of fitted() (cumulative: P(y = 1)) at x = 1.2, 1.5, 1.8
# over the poly() fit's.
#   ILLSE_ARM=base|lane Rscript dev/illse-glm.R [seeds] [families]
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
args <- commandArgs(TRUE)
seeds <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:3
fams <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
  c("bernoulli", "poisson", "cumulative")
tag <- function(w) {
  if (!length(w)) return("-")
  k <- c(illcond = "ill-conditioned", short = "stopped short",
         lost = "Standard errors are not available",
         grad = "gradient", flat = "flat along", conv = "false conv",
         sep = "separat")
  paste(names(k)[vapply(k, function(s) any(grepl(s, w)), NA)],
        collapse = ",")
}
mpm <- function(A) mpfrArray(A, mp_prec, dim = dim(A))
# B with X_raw = X_poly B over the coefficient columns, exactly
basis_change <- function(Xr, Xp) {
  Pm <- mpm(Xp)
  mp_solve(crossprod(Pm)) %*% crossprod(Pm, mpm(Xr))
}
exact_glm <- function(fam, X, y, p) {
  b <- p[names(p) == "beta"]
  Xm <- mpm(X)
  eta <- Xm %*% mp(b)
  if (fam == "bernoulli") {
    pr <- 1 / (1 + exp(-eta))
    w <- pr * (1 - pr)
  } else {
    w <- exp(eta)
  }
  mp_solve(crossprod(Xm * w[, 1], Xm))
}
exact_cum <- function(X, y, p) {
  k <- as.integer(y)
  nb <- sum(names(p) == "beta")
  Xm <- mpm(X)
  nll <- function(th) {
    b <- th[seq_len(nb)]
    raw <- th[-seq_len(nb)]
    tau <- c(raw[1], raw[1] + cumsum(exp(raw[-1])))
    eta <- Xm %*% b
    cdf <- function(z) 1 / (1 + exp(-z))
    K <- length(tau) + 1
    hi <- mp(rep(1, length(k)))
    lo <- mp(rep(0, length(k)))
    up <- which(k < K)
    dn <- which(k > 1)
    hi[up] <- cdf(tau[k[up]] - eta[up])
    lo[dn] <- cdf(tau[k[dn] - 1] - eta[dn])
    -sum(log(hi - lo))
  }
  mp_solve(mp_hess(nll, p))
}
xr <- c(1.2, 1.5, 1.8)
for (fam in fams) for (deg in 6:8) for (re in c(FALSE, TRUE)) {
  for (seed in seeds) {
    ds <- des_poly(deg, re, fam)(seed)
    r <- cap(frm(ds$formula, data = ds$data, family = ds$family))
    if (inherits(r$v, "illse_err")) {
      cat(sprintf("%-26s s%d ERROR %s\n", ds$lab, seed, r$v)); next
    }
    f <- r$v
    g <- suppressWarnings(frm(ds$poly_formula, data = ds$data,
                              family = ds$family))
    p <- f$opt$par
    pb <- which(names(p) == "beta")
    X <- as.matrix(f$frame$linpreds[[1]]$X)
    Xp <- as.matrix(g$frame$linpreds[[1]]$X)
    V <- suppressWarnings(vcov(f, full = TRUE))
    Vb <- V[pb, pb, drop = FALSE]
    # the poly() fit's coefficient block, carried to raw coefficients
    Vg <- suppressWarnings(vcov(g, full = TRUE))
    pg <- which(names(g$opt$par) == "beta")
    B <- basis_change(X, Xp)
    Bi <- mp_solve(B)
    Vp <- Bi %*% mpm(Vg[pg, pg, drop = FALSE]) %*% t(Bi)
    nan3 <- c(cor = NaN, rel = NaN, diag = NaN)
    ep <- if (all(is.finite(Vb))) mp_cov_err(Vb, Vp) else nan3
    ee <- nan3
    if (!re && all(is.finite(V))) {
      Ve <- if (fam == "cumulative") exact_cum(X, ds$data$y, p) else
        exact_glm(fam, X, ds$data$y, p)
      ee <- mp_cov_err(if (fam == "cumulative") V else Vb, Ve)
    }
    nd <- ds$data[rep(1, length(xr)), , drop = FALSE]
    nd$x <- xr
    q <- function(m) {
      out <- cap(fitted(m, newdata = nd, re_formula = NA))
      if (inherits(out$v, "illse_err")) return(list(se = rep(NaN, 3),
                                                    w = out$w))
      v <- out$v
      se <- if (length(dim(v)) == 3L) v[, "Est.Error", 1] else
        v[, "Est.Error"]
      list(se = se, w = out$w)
    }
    qf <- q(f)
    qg <- q(g)
    cat(sprintf(paste0("%-26s s%d code %d | warn %-14s | lost %-10s | ",
                       "steps %s | ll gap %9.2e | vs exact cor %.1e rel ",
                       "%.1e diag %.1e | vs poly cor %.1e diag %.1e | ",
                       "fitted SE/poly %s %s\n"),
                ds$lab, seed, f$opt$convergence, tag(r$w), lost_of(f),
                if (is.null(f$opt$polish)) 0 else f$opt$polish[["steps"]],
                as.numeric(logLik(g)) - as.numeric(logLik(f)),
                ee[["cor"]], ee[["rel"]], ee[["diag"]], ep[["cor"]],
                ep[["diag"]],
                paste(sprintf("%.4f", qf$se / qg$se), collapse = " "),
                tag(qf$w)))
  }
}
