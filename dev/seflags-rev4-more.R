# Reviewer of lane seflags, check of punch 2b: more variance routes on
# ill-conditioned fits, each against the same call on the poly(x, 8)
# fit (same likelihood, well conditioned).
#  1. hypothesis() of the polynomial at x = 1.5, 2 and 2.5 (precomputed
#     columns x1..x8), no random effects, seeds 1 to 12 (hyp_par_cov's
#     ML branch, whose `ill` needs a non-NULL se_null).
#  2. REML with (1 | g): fitted() and vcov() (beta is integrated out, so
#     the check never sees the direction), lane against base.
#  3. bernoulli and cumulative with the raw polynomial: fitted() on the
#     response scale (the cumulative one through fit_fd_se()).
#  4. hypothesis(vcov = vcov_cluster()) and a user's own a' vcov() a.
#   Rscript dev/seflags-rev4-more.R <lane|base|op4>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-common.R")
rat <- function(a, b) {
  if (inherits(a, "rev_err") || inherits(b, "rev_err")) {
    return(paste("ERR", substr(c(a, b)[1], 1, 90)))
  }
  paste(format(as.numeric(a) / as.numeric(b), digits = 4), collapse = " ")
}
mkd <- function(s, re_sd = 0) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$eta <- sin(3 * d$x) + stats::rnorm(20, 0, re_sd)[d$g]
  d$y <- d$eta + stats::rnorm(200, 0, 0.3)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  d
}
rhs <- paste0("x", 1:8, collapse = " + ")
hy <- function(x0) {
  paste0("Intercept + ", paste0(sprintf("%.17g", x0^(1:8)), " * x", 1:8,
                                collapse = " + "), " = 0")
}
nd <- data.frame(x = c(1.5, 2, 2.5), g = factor(1))
for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
cat("\n## 1. hypothesis() of the polynomial, no random effects\n")
for (s in 1:12) {
  d <- mkd(s)
  r <- cap(frm(stats::as.formula(paste("y ~", rhs)), data = d))
  f <- r$v
  if (f$opt$convergence != 0 || !any(grepl("ill-conditioned", r$w))) {
    cat(sprintf("seed %2d | code %d | ill %s | skipped\n", s,
                f$opt$convergence, any(grepl("ill-conditioned", r$w))))
    next
  }
  g <- suppressWarnings(frm(y ~ poly(x, 8), data = d))
  ndp <- nd
  ref <- suppressWarnings(fitted(g, newdata = ndp)[, "Est.Error"])
  hs <- vapply(nd$x, function(x0) {
    h <- cap(hypothesis(f, hy(x0)))
    if (inherits(h$v, "rev_err")) NA_real_ else h$v$hypothesis$Est.Error
  }, 0)
  fi <- cap(fitted(f, newdata = nd))$v[, "Est.Error"]
  sdr <- suppressWarnings(ns$sdr_of(f))
  cat(sprintf("seed %2d | se_null %s | hypothesis/ref %s | fitted/ref %s\n", s,
              if (is.null(sdr$se_null)) "NULL" else paste(dim(sdr$se_null), collapse = "x"),
              rat(hs, ref), rat(fi, ref)))
}
cat("\n## 2. REML with (1 | g)\n")
for (s in c(3, 8, 10)) {
  d <- mkd(s, 0.3)
  r <- cap(frm(stats::as.formula(paste("y ~", rhs, "+ (1 | g)")), data = d,
               REML = TRUE))
  f <- r$v
  if (inherits(f, "rev_err")) {
    cat(sprintf("seed %2d | REML fit failed: %s
", s, substr(f, 1, 90)))
    next
  }
  g <- suppressWarnings(frm(y ~ poly(x, 8) + (1 | g), data = d, REML = TRUE))
  a <- cap(fitted(f, newdata = nd, re_formula = NA))
  b <- suppressWarnings(fitted(g, newdata = nd, re_formula = NA))
  cat(sprintf("seed %2d | code %d (ref %d) | warnings %d%s | fitted/ref %s (warn %d)\n",
              s, f$opt$convergence, g$opt$convergence, length(r$w),
              if (length(r$w)) paste0(" (", substr(r$w[1], 1, 50), ")") else "",
              rat(a$v[, "Est.Error"], b[, "Est.Error"]), length(a$w)))
}
cat("\n## 3. bernoulli and cumulative\n")
for (s in c(2, 5, 8)) {
  d <- mkd(s)
  d$yb <- as.integer(d$y > stats::median(d$y))
  d$yo <- cut(d$y, stats::quantile(d$y, c(0, 0.3, 0.7, 1)),
              include.lowest = TRUE, labels = FALSE)
  for (fam in c("bernoulli", "cumulative")) {
    lhs <- if (fam == "bernoulli") "yb" else "yo"
    fm <- if (fam == "bernoulli") bernoulli() else cumulative()
    r <- cap(frm(stats::as.formula(paste(lhs, "~", rhs)), data = d,
                 family = fm))
    f <- r$v
    g <- suppressWarnings(frm(stats::as.formula(paste(lhs, "~ poly(x, 8)")),
                              data = d, family = fm))
    a <- cap(fitted(f, newdata = nd))
    b <- cap(fitted(g, newdata = nd))
    pick <- function(z) {
      if (inherits(z, "rev_err")) return(z)
      if (is.matrix(z)) z[, "Est.Error"] else as.vector(z[, "Est.Error", ])
    }
    cat(sprintf("seed %d %-10s | code %d (ref %d) | ill %s | other warnings %d | fitted/ref %s (warn %d)\n",
                s, fam, f$opt$convergence, g$opt$convergence,
                any(grepl("ill-conditioned", r$w)),
                sum(!grepl("ill-conditioned", r$w)),
                rat(pick(a$v), pick(b$v)), length(a$w)))
  }
}
cat("\n## 4. vcov_cluster() and a user's own a' vcov() a\n")
for (s in c(2, 5, 9)) {
  d <- mkd(s)
  d$cl <- factor(rep(1:40, each = 5))
  f <- suppressWarnings(frm(stats::as.formula(paste("y ~", rhs)), data = d))
  g <- suppressWarnings(frm(y ~ poly(x, 8), data = d))
  if (f$opt$convergence != 0) { cat("seed", s, "code", f$opt$convergence, "\n"); next }
  a <- cap(hypothesis(f, hy(2), vcov = vcov_cluster(f, cluster = d$cl)))
  b <- cap(hypothesis(g, "Intercept = 0"))
  pa <- stats::model.matrix(~ poly(x, 8), d)  # unused, keeps the bases
  A0 <- c(1, 2^(1:8))
  V <- stats::vcov(f)
  own <- sqrt(drop(t(A0) %*% V %*% A0))
  fi <- suppressWarnings(fitted(f, newdata = nd[2, ]))[, "Est.Error"]
  cat(sprintf("seed %d | hypothesis(vcov = vcov_cluster()) at x = 2: %s (warn %d) | a' vcov() a / fitted(): %s\n",
              s, if (inherits(a$v, "rev_err")) substr(a$v, 1, 60) else
                format(a$v$hypothesis$Est.Error, digits = 4),
              length(a$w), rat(own, fi)))
}
