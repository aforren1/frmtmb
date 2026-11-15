# Lane nlstop2, item 5053: on the stops of dev/nlstop2-ridgeprior.R moved
# along the ridge of c0 + exp(a)^k (psd 100 k 5, psd 1000 k 20, 50, 80),
# what the quadratic model predicts against what Newton steps on the
# stop's Hessian (chord iterations, the Hessian held) actually gain.
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-ridgedbg.R
source("dev/nlstop2-common.R")
ns <- asNamespace("frmtmb")
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
for (cs in list(c(100, 5), c(1000, 20), c(1000, 50), c(1000, 80))) {
  psd <- cs[1]; kt <- cs[2]
  pr <- set_prior(sprintf("normal(0.5, %s)", psd), nlpar = "k")
  f <- suppressWarnings(suppressMessages(frm(form, data = dn, prior = pr)))
  nm <- ns$outer_par_names(f)
  ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
  it <- match("theta_1", nm)
  opt <- function(par, fn, gr, lower, upper, control) {
    o <- stats::nlminb(par, fn, gr, control = list(iter.max = 1000,
                                                   eval.max = 2000))
    p <- o$par
    s <- p[ik] / kt
    p[ia] <- p[ia] * s; p[ik] <- kt; p[it] <- p[it] + log(abs(s))
    list(par = p, objective = fn(p), convergence = 0L,
         message = "relative convergence (4)")
  }
  r <- tt(frm(form, data = dn, prior = pr,
              control = frmtmb_control(optimizer = opt, restarts = 0)))
  obj <- r$obj; p <- r$opt$par
  g <- drop(obj$gr(p))
  # the Hessian grad_headroom() reads
  u1 <- is.null(r$par_units) || all(r$par_units == 1)
  H <- ns$fit_outer_hessian_at(r, p, kept_only = !u1)
  if (is.null(H)) {
    H <- ns$fd_hessian(p, function(q) drop(obj$gr(q)),
                       eps = ns$outer_fd_steps(r$frame, names(p)))$H
  }
  H <- (H + t(H)) / 2
  cat("headroom (grad_headroom):", ns$grad_headroom(r, g, rep(FALSE, length(p))),
      "\n")
  # the coordinates grad_headroom() takes (the design's, refined)
  h <- ns$headroom_reparam_hessian(r, H, seq_along(p))
  if (is.null(h)) h <- list(H = H, M = diag(length(p)))
  ch <- chol(h$H)
  M <- h$M
  pred <- 0.5 * sum(backsolve(ch, drop(crossprod(M, g)),
                              transpose = TRUE)^2)
  short <- r$opt$objective - f$opt$objective
  cat(sprintf("psd %g k %g: short %.4g gmax %.3g pred %.4g warn %s\n",
              psd, kt, short, max(abs(g)), pred, wsum(r)))
  # the same stop as a start of the default optimizer (nlminb, restarts)
  est <- r$estimates
  r2 <- tt(frm(form, data = dn, prior = pr,
               start = est[intersect(c("beta", "betad", "theta"), names(est))]))
  cat(sprintf("  default optimizer from the stop: short %.4g code %s warn %s\n",
              r2$opt$objective - f$opt$objective, r2$opt$convergence,
              wsum(r2)))
  r3 <- tt(frm(form, data = dn, prior = pr,
               control = frmtmb_control(restarts = 0),
               start = est[intersect(c("beta", "betad", "theta"), names(est))]))
  cat(sprintf("  nlminb from the stop, restarts = 0: short %.4g code %s warn %s\n",
              r3$opt$objective - f$opt$objective, r3$opt$convergence,
              wsum(r3)))
  q <- p; f0 <- obj$fn(p)
  for (i in 1:10) {
    gq <- drop(crossprod(M, drop(obj$gr(q))))
    q <- q - drop(M %*% backsolve(ch, backsolve(ch, gq, transpose = TRUE)))
    cat(sprintf("  chord %d: gain %.4g\n", i, f0 - obj$fn(q)))
  }
  # Newton with the Hessian retaken at each step
  q <- p
  for (i in 1:4) {
    gq <- drop(obj$gr(q))
    Hq <- stats::optimHess(q, function(z) obj$fn(z), function(z) obj$gr(z))
    st <- tryCatch(solve((Hq + t(Hq)) / 2, gq), error = function(e) NULL)
    if (is.null(st)) break
    q <- q - st
    cat(sprintf("  newton %d: gain %.4g\n", i, f0 - obj$fn(q)))
  }
}
cat("DONE\n")
