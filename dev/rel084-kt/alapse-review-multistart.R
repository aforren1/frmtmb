# rel084-kt copy of dev/reviews/alapse-review-multistart.R: sources dev/rel084-kt/alapse-review-common.R (rellib-r23). No other change.
# Review of lane alapse, check 1: on designs D (sigma ~ z, quantile
# estimated constant) and E (quantile ~ z) of dev/alapse-rq.R, is the
# proved fit the highest maximum? Each fit is refitted with the quantile
# fixed on a grid (quantile = p0 for E's intercept start is not possible
# through frm, so E uses alternating rq()/optim() from grid starts).
#   Rscript dev/reviews/alapse-review-multistart.R lane <seeds>
a <- commandArgs(TRUE)
lib_arg <- a[1]
seeds <- eval(parse(text = a[2]))
source("dev/rel084-kt/alapse-review-common.R")

rq_b <- function(y, X, w, p) {
  if (length(unique(p)) == 1L) {
    return(suppressWarnings(quantreg::rq.wfit(X, y, tau = p[1], weights = w,
                                              method = "br"))$coefficients)
  }
  cc <- colSums(w * (2 * p - 1) * X)
  M <- 1e4 * (1 + max(abs(y)))
  suppressWarnings(quantreg::rq.wfit(rbind(X, cc), c(y, M), tau = 0.5,
                                     weights = c(w, 1),
                                     method = "br"))$coefficients
}
alt_max <- function(y, X, Xs, Xq, bs, bq) {
  ll_of <- function(b, bs, bq) {
    ald_ll(y, drop(X %*% b), exp(drop(Xs %*% bs)), plogis(drop(Xq %*% bq)))
  }
  cur <- -Inf
  for (it in 1:300) {
    b <- rq_b(y, X, 1 / exp(drop(Xs %*% bs)), plogis(drop(Xq %*% bq)))
    ks <- length(bs)
    o <- optim(c(bs, bq), function(t) -ll_of(b, t[seq_len(ks)],
                                             t[-seq_len(ks)]),
               method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
    bs <- o$par[seq_len(ks)]
    bq <- o$par[-seq_len(ks)]
    new <- ll_of(b, bs, bq)
    if (abs(new - cur) < 1e-11) { cur <- new; break }
    cur <- new
  }
  list(ll = cur, bq = bq)
}

for (des in c("D", "E")) {
  res <- NULL
  for (s in seeds) {
    if (des == "D") {
      set.seed(4000 + s)
      d <- data.frame(x = rnorm(200), z = rnorm(200))
      d$y <- 1 + 0.5 * d$x + exp(0.4 * d$z) * rnorm(200)
      fo <- bf(y ~ x, sigma ~ z)
      Xs <- cbind(1, d$z); Xq <- cbind(rep(1, 200))
    } else {
      set.seed(5000 + s)
      d <- data.frame(x = rnorm(200), z = rnorm(200))
      d$y <- 1 + 0.5 * d$x + rnorm(200)
      fo <- bf(y ~ x, quantile ~ z)
      Xs <- cbind(rep(1, 200)); Xq <- cbind(1, d$z)
    }
    X <- cbind(1, d$x)
    r <- cap(frm(fo, family = asym_laplace(), data = d))
    f <- r$v
    ll <- as.numeric(logLik(f))
    best <- -Inf
    bbq <- NA
    for (q0 in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
      bq0 <- c(qlogis(q0), rep(0, ncol(Xq) - 1))
      am <- alt_max(d$y, X, Xs, Xq, rep(0, ncol(Xs)), bq0)
      if (am$ll > best) { best <- am$ll; bbq <- am$bq[1] }
    }
    conv <- any(r$wc %in% c("frmtmb_conv_code", "frmtmb_conv_gradient"))
    cat(sprintf("%s %3d code %d cert %-5s conv-warn %-5s ll %.6f multistart %.6f gap %.3e q-int %.4f best q-int %.4f\n",
                des, s, f$opt$convergence, isTRUE(f$opt$ald_cert$ok), conv, ll,
                best, best - ll,
                f$opt$par[grep("quantile", ns$outer_par_names(f))[1]], bbq))
    res <- rbind(res, data.frame(gap = best - ll,
                                 cert = isTRUE(f$opt$ald_cert$ok),
                                 conv = conv))
  }
  cat(sprintf("%s summary: n %d | proved %d | proved, gap > 1e-6, silent: %d | max gap %.3g\n",
              des, nrow(res), sum(res$cert),
              sum(res$cert & res$gap > 1e-6 & !res$conv), max(res$gap)))
}
cat("DONE\n")
