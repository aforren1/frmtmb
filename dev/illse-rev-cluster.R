# Reviewer of lane illse, attack 3: vcov_cluster() on a degree-8 raw
# polynomial with (1 | g): the sandwich B S'S B formed in double from
# (a) the lane's vcov(full = TRUE), (b) the exact bread rounded to
# double, (c) through the lane's factor, root (root' S'S root) root',
# each against the exact sandwich (256 bits, same scores S).
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
deg <- 8; xn <- paste0("x", 1:deg)
for (seed in 101:103) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(200, 0, 0.3)
  for (k in 1:deg) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
  f <- suppressWarnings(frm(fo, data = d))
  p <- f$opt$par; pn <- names(p)
  pb <- which(pn == "beta"); pu <- which(pn == "betad"); pw <- which(pn == "theta")
  X <- as.matrix(f$frame$linpreds[[1]]$X)
  Xm <- mM(X); ym <- m(d$y)
  XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
  gi <- split(seq_len(nrow(d)), d$g)
  S <- lapply(gi, function(ix) list(s = colSums(Xm[ix, , drop = FALSE]), sy = sum(ym[ix]), n = length(ix)))
  nll <- function(th) {
    b <- th[pb]; a2 <- exp(2 * th[pu]); c2 <- exp(2 * th[pw])
    bm <- mpfr2array(b, dim = c(length(pb), 1L))
    out <- (yy - 2 * sum(Xy * b) + sum(bm * (XX %*% bm))) / a2
    for (st in S) {
      sr <- st$sy - sum(st$s * b)
      out <- out + (st$n - 1) * log(a2) + log(a2 + st$n * c2) - c2 * sr^2 / (a2 * (a2 + st$n * c2))
    }
    out / 2
  }
  Ve <- mp_inv(mp_hessian(nll, p))
  Sc <- ns$cluster_scores_at(f, ns$resolve_cluster(f, ~ g))
  SS <- crossprod(Sc)
  Vce <- Ve %*% crossprod(mM(Sc)) %*% Ve
  V <- suppressWarnings(vcov(f, full = TRUE))
  root <- ns$sdr_of(f)$cov_fixed_root
  Va <- V %*% SS %*% V
  Vb <- num(Ve) %*% SS %*% num(Ve)
  K <- Sc %*% root
  Vc <- root %*% crossprod(K) %*% t(root)
  vc <- suppressWarnings(vcov_cluster(f, ~ g, full = TRUE))
  A <- rbind(c(1, 1.5^(1:deg), 0, 0))
  ex <- sqrt(mp_rowvar(A, Vce))
  sev <- function(M) sqrt(max(drop(A %*% M %*% t(A)), 0)) / ex
  sef <- if (TRUE) sqrt(sum((A %*% root %*% t(K))^2)) / ex
  cat(sprintf(paste0("s%d | entries cor-scale error: vcov_cluster %.2e | (a) lane V whole %.2e | ",
                     "(b) exact bread rounded %.2e | (c) through factor %.2e | ",
                     "SE of fitted at 1.5 / exact: vcov_cluster %.4g, (c) whole %.4g, (c) as factor %.6f\n"),
              seed, cov_err(vc, Vce)[["cor"]], cov_err(Va, Vce)[["cor"]], cov_err(Vb, Vce)[["cor"]],
              cov_err(Vc, Vce)[["cor"]], sev(unclass(vc)), sev(Vc), sef))
}
