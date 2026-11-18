# Lane alapse, backlog 4982: each asym_laplace() fit without random
# effects against quantreg::rq() at the fit's own quantile, and against
# the exact joint maximum.
#   Rscript dev/alapse-rq.R <base|lane> [seeds, default 1:25]
# Designs (n, formula, noise):
#   A  160, y ~ x, ALD noise (p 0.3)                  quantile estimated
#   B  100, y ~ x, gaussian noise                     quantile estimated
#   C  120, y ~ x1 + x2, t3 noise, quantile = 0.25    quantile fixed
#   D  200, y ~ x, sigma ~ z, gaussian noise scaled by exp(0.4 z)
#   E  200, y ~ x, quantile ~ z, gaussian noise
#   F  200, zero_inflated_asym_laplace, y ~ x, 20% zeros
#   G  300, y ~ x + fac (3 levels), quantile estimated, ALD noise (p 0.7)
# plus the 4982 cases of dev/selimit2-asym.R (dbg14, gauss5, fixq, zi).
# Per fit: frm's code, warning classes, logLik; `gap_b` = the ALD
# log-likelihood at rq()'s coefficients (weights 1 / sigma_i, the fit's
# own quantiles) minus frm's, the fit's sigma and quantile held; `gap_j`
# = the joint maximum (the quantile process for a constant sigma and
# quantile, else alternating rq() and optim() from frm's estimates)
# minus frm's logLik; `db` = max |frm coefficient - rq coefficient|.
a <- commandArgs(TRUE)
lib_arg <- a[1]
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:25
source("dev/alapse-common.R")
suppressMessages(library(quantreg))

# b-step: the weighted check-function minimum at row quantiles p, as a
# median regression with one pseudo-row carrying the linear part
rq_b <- function(y, X, w, p) {
  if (length(unique(p)) == 1L) {
    fit <- quantreg::rq.wfit(X, y, tau = p[1], weights = w, method = "br")
    return(fit$coefficients)
  }
  cc <- colSums(w * (2 * p - 1) * X)
  M <- 1e4 * (1 + max(abs(y)))
  Xa <- rbind(X, cc)
  ya <- c(y, M)
  fit <- quantreg::rq.wfit(Xa, ya, tau = 0.5, weights = c(w, 1),
                           method = "br")
  fit$coefficients
}

# the exact joint maximum for a constant sigma and quantile: the
# quantile process holds every distinct rq() solution; on each interval
# the profile in p is smooth
exact_const <- function(y, X, p_fixed = NULL) {
  n <- length(y)
  if (!is.null(p_fixed)) {
    b <- rq_b(y, X, rep(1, n), rep(p_fixed, n))
    r <- y - drop(X %*% b)
    s <- mean(r * (p_fixed - (r < 0)))
    return(list(ll = ald_ll(y, drop(X %*% b), s, p_fixed), b = b, p = p_fixed,
                s = s))
  }
  pr <- quantreg::rq(y ~ X - 1, tau = -1)
  sol <- pr$sol
  taus <- sol[1, ]
  B <- sol[-(1:3), , drop = FALSE]
  best <- list(ll = -Inf)
  for (j in seq_len(ncol(B))) {
    lo <- taus[j]
    hi <- if (j < ncol(B)) taus[j + 1] else 1
    lo <- max(lo, 1e-9); hi <- min(hi, 1 - 1e-9)
    if (hi <= lo) next
    r <- y - drop(X %*% B[, j])
    A <- sum(pmax(r, 0)); Bn <- sum(pmax(-r, 0))
    prof <- function(p) n * log(p) + n * log1p(-p) - n * log((p * A +
                                    (1 - p) * Bn) / n) - n
    o <- optimize(prof, c(lo, hi), maximum = TRUE, tol = 1e-12)
    for (pp in c(o$maximum, lo, hi)) {
      v <- prof(pp)
      if (is.finite(v) && v > best$ll) {
        best <- list(ll = v, b = B[, j], p = pp,
                     s = (pp * A + (1 - pp) * Bn) / n)
      }
    }
  }
  best
}

# alternating maximum from a start: rq() step, then optim() on the
# sigma and quantile coefficients
alt_max <- function(y, X, Xs, Xq, b, bs, bq, p_fixed = NULL, w = 1) {
  ll_of <- function(b, bs, bq) {
    s <- exp(drop(Xs %*% bs))
    p <- if (is.null(p_fixed)) plogis(drop(Xq %*% bq)) else
      rep(p_fixed, length(y))
    ald_ll(y, drop(X %*% b), s, p, w)
  }
  cur <- ll_of(b, bs, bq)
  for (it in 1:200) {
    s <- exp(drop(Xs %*% bs))
    p <- if (is.null(p_fixed)) plogis(drop(Xq %*% bq)) else
      rep(p_fixed, length(y))
    b <- rq_b(y, X, w / s, p)
    th <- c(bs, if (is.null(p_fixed)) bq)
    ks <- length(bs)
    o <- optim(th, function(t) -ll_of(b, t[seq_len(ks)],
                                     if (is.null(p_fixed)) t[-seq_len(ks)]),
               method = "BFGS", control = list(reltol = 1e-14,
                                               maxit = 1000))
    bs <- o$par[seq_len(ks)]
    if (is.null(p_fixed)) bq <- o$par[-seq_len(ks)]
    new <- ll_of(b, bs, bq)
    if (abs(new - cur) < 1e-12) { cur <- new; break }
    cur <- new
  }
  list(ll = cur, b = b, bs = bs, bq = bq, it = it)
}

# zi data's ALD rows
one_case <- function(id, seed, d, fo, fam, p_fixed = NULL, zi = FALSE,
                     sform = ~1, qform = ~1) {
  r <- cap(frm(fo, family = fam, data = d))
  if (inherits(r$v, "error")) {
    cat(sprintf("%s %s ERROR %s\n", id, seed, conditionMessage(r$v)))
    return(invisible())
  }
  f <- r$v
  nm <- ns$outer_par_names(f)
  par <- setNames(f$opt$par, nm)
  dd <- if (zi) d[d$y != 0, , drop = FALSE] else d
  X <- model.matrix(fo_mu(fo), dd)
  Xs <- model.matrix(sform, dd)
  Xq <- model.matrix(qform, dd)
  b <- par[colnames(X)]
  bs <- par[paste0("sigma_", colnames(Xs))]
  bq <- if (is.null(p_fixed)) par[paste0("quantile_", colnames(Xq))]
  s <- exp(drop(Xs %*% bs))
  p <- if (is.null(p_fixed)) plogis(drop(Xq %*% bq)) else
    rep(p_fixed, nrow(dd))
  y <- dd$y
  ll_frm_ald <- ald_ll(y, drop(X %*% b), s, p)
  brq <- rq_b(y, X, 1 / s, p)
  gap_b <- ald_ll(y, drop(X %*% brq), s, p) - ll_frm_ald
  ex <- if (identical(sform, ~1) && identical(qform, ~1)) {
    exact_const(y, X, p_fixed)
  } else {
    alt_max(y, X, Xs, Xq, b, bs, bq, p_fixed)
  }
  gap_j <- ex$ll - ll_frm_ald
  # logLik of frm against the ALD rows' own sum: the zi part is apart
  ll <- as.numeric(logLik(f))
  cls <- unique(c(r$wc, r$mc))
  cls <- sub("^frmtmb_", "", cls)
  cat(sprintf(paste0("%-6s %3s code %d %-26s ll %.8f | gap_b %9.2e | ",
                     "gap_j %9.2e | db %8.2e | dp %8.2e | %s\n"),
              id, seed, f$opt$convergence,
              substr(f$opt$message, 1, 26), ll, gap_b, gap_j,
              max(abs(b - brq)),
              if (is.null(p_fixed) && identical(qform, ~1)) {
                if (is.null(ex$p)) NA_real_ else abs(p[1] - ex$p)
              } else NA_real_,
              paste(cls, collapse = ",")))
  invisible()
}
fo_mu <- function(fo) {
  if (inherits(fo, "brmsformula") || inherits(fo, "frmtmb_bf") ||
        is.list(fo)) {
    g <- fo$formula %||% fo[[1]]
  } else g <- fo
  stats::reformulate(attr(stats::terms(g), "term.labels"), intercept = TRUE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# the 4982 cases
set.seed(43); n <- 160; x <- rnorm(n)
y <- brms::rskew_normal(n, 1 + 0.5 * x, 1, 4)
y <- brms::rexgaussian(n, 2 + 0.5 * x, 0.5, 1)
y <- brms::rasym_laplace(n, 1 + 0.5 * x, 1, 0.3)
one_case("dbg14", 43, data.frame(y = y, x = x), y ~ x, asym_laplace())
set.seed(5)
d <- data.frame(x = stats::rnorm(100))
d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
one_case("gauss5", 5, d, y ~ x, asym_laplace())
one_case("fixq", 5, d, bf(y ~ x, quantile = 0.5), asym_laplace(),
         p_fixed = 0.5)
set.seed(9)
dz <- data.frame(x = stats::rnorm(200))
dz$y <- ifelse(stats::runif(200) < 0.2, 0,
               brms::rasym_laplace(200, 1 + 0.5 * dz$x, 1, 0.5))
one_case("zi", 9, dz, y ~ x, zero_inflated_asym_laplace(), zi = TRUE)

for (s in seeds) {
  set.seed(1000 + s)
  n <- 160; x <- rnorm(n)
  d <- data.frame(x = x, y = brms::rasym_laplace(n, 1 + 0.5 * x, 1, 0.3))
  one_case("A", s, d, y ~ x, asym_laplace())
  set.seed(2000 + s)
  d <- data.frame(x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  one_case("B", s, d, y ~ x, asym_laplace())
  set.seed(3000 + s)
  d <- data.frame(x1 = rnorm(120), x2 = runif(120))
  d$y <- 2 - d$x1 + 3 * d$x2 + rt(120, 3)
  one_case("C", s, d, bf(y ~ x1 + x2, quantile = 0.25), asym_laplace(),
           p_fixed = 0.25)
  set.seed(4000 + s)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- 1 + 0.5 * d$x + exp(0.4 * d$z) * rnorm(200)
  one_case("D", s, d, bf(y ~ x, sigma ~ z), asym_laplace(), sform = ~z)
  set.seed(5000 + s)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- 1 + 0.5 * d$x + rnorm(200)
  one_case("E", s, d, bf(y ~ x, quantile ~ z), asym_laplace(), qform = ~z)
  set.seed(6000 + s)
  d <- data.frame(x = rnorm(200))
  d$y <- ifelse(runif(200) < 0.2, 0,
                brms::rasym_laplace(200, 1 + 0.5 * d$x, 1, 0.5))
  one_case("F", s, d, y ~ x, zero_inflated_asym_laplace(), zi = TRUE)
  set.seed(7000 + s)
  d <- data.frame(x = rnorm(300), fac = factor(sample(c("a", "b", "c"), 300,
                                                      TRUE)))
  d$y <- brms::rasym_laplace(300, 1 + 0.5 * d$x + (d$fac == "b") -
                               0.5 * (d$fac == "c"), 1, 0.7)
  one_case("G", s, d, y ~ x + fac, asym_laplace())
}
cat("DONE\n")
