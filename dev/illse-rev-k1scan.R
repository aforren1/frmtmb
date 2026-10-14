# Reviewer of lane illse, attack 4: where the first pass in the design's
# coordinates lands relative to se_reparam_refine (1e4), and what the
# refinement changes. For each fit that takes the route: the first-pass
# kappa (unit diagonal), whether the refinement ran, and the refined
# against the unrefined covariance on the correlation scale.
#   Rscript dev/illse-rev-k1scan.R <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
seeds <- eval(parse(text = commandArgs(TRUE)[1]))
mk <- function(seed, deg) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
  u <- stats::rnorm(20, 0, 0.5)[d$g]
  e <- sin(3 * d$x)
  d$yg <- e + u + stats::rnorm(n, 0, 0.3)
  d$yb <- stats::rbinom(n, 1, stats::plogis(2 * e + u))
  d$yp <- stats::rpois(n, exp(0.5 + e + u))
  d$yz <- ifelse(stats::runif(n) < stats::plogis(-1 + e), 0, stats::rpois(n, exp(1 + e)))
  d$yo <- factor(cut(2 * e + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE),
                 ordered = TRUE)
  d$ys <- e + stats::rnorm(n, 0, exp(-1 + 0.5 * e))
  d$yt <- e + 0.3 * stats::rt(n, 5)
  for (k in seq_len(deg)) d[[paste0("x", k)]] <- d$x^k
  d
}
rhs <- function(deg) paste(paste0("x", seq_len(deg)), collapse = " + ")
specs <- function(deg) {
  r <- rhs(deg)
  list(
    gauss_re = list(stats::as.formula(paste("yg ~", r, "+ (1 | g)")), gaussian()),
    bern = list(stats::as.formula(paste("yb ~", r)), bernoulli()),
    bern_re = list(stats::as.formula(paste("yb ~", r, "+ (1 | g)")), bernoulli()),
    pois_re = list(stats::as.formula(paste("yp ~", r, "+ (1 | g)")), poisson()),
    cum = list(stats::as.formula(paste("yo ~", r)), cumulative()),
    cum_re = list(stats::as.formula(paste("yo ~", r, "+ (1 | g)")), cumulative()),
    zip = list(bf(stats::as.formula(paste("yz ~", r)), zi ~ x), zero_inflated_poisson()),
    sigma = list(bf(stats::as.formula(paste("ys ~", r)),
                    stats::as.formula(paste("sigma ~", r))), gaussian()),
    student = list(stats::as.formula(paste("yt ~", r)), student()))
}
for (seed in seeds) for (deg in c(5, 7)) {
  d <- mk(seed, deg)
  for (nm in names(specs(deg))) {
    sp <- specs(deg)[[nm]]
    r <- cap(frm(sp[[1]], data = d, family = sp[[2]]))
    f <- r$v
    if (inherits(f, "rev_err")) {
      cat(sprintf("%-9s deg %d s%d ERROR %s\n", nm, deg, seed, substr(f, 1, 60))); next
    }
    root <- ns$sdr_of(f)$cov_fixed_root
    p <- f$opt$par
    M <- tryCatch(ns$reparam_basis(f$frame, names(p)), error = function(e) NULL)
    k1 <- NA; refined <- NA; dref <- NA
    if (!is.null(root) && !is.null(M)) {
      h1 <- ns$reparam_hessian(f$obj, p, M)
      D1 <- sqrt(abs(diag(h1$H)))
      e1 <- eigen(h1$H / outer(D1, D1), symmetric = TRUE, only.values = TRUE)$values
      k1 <- max(e1) / min(e1)
      h2 <- ns$reparam_hessian_refined(f$obj, p, M)
      refined <- !identical(h2$M, M)
      V1 <- tryCatch(M %*% solve(h1$H) %*% t(M), error = function(e) NULL)
      V2 <- h2$M %*% solve(h2$H) %*% t(h2$M)
      dg <- sqrt(diag(V2))
      dref <- if (is.null(V1)) NA else max(abs(V1 - V2) / outer(dg, dg))
      f$obj$fn(p)
    }
    cat(sprintf(paste0("%-9s deg %d s%d | code %d | route %s lossy %s | lost %d | ",
                       "first-pass kappa %.2e | refined %s | refined vs unrefined %.1e | warn %s\n"),
                nm, deg, seed, f$opt$convergence, !is.null(root),
                isTRUE(attr(root, "lossy")), length(ns$sdr_of(f)$se_lost),
                k1, refined, dref, wtag(r$w)))
  }
}
