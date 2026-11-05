# Fix round 1, item 1: spellings beside the warm start. My own dense marginal, no frmtmb code, each with
# many starts: two exact blocks, by = f, by = w, weights (with the weighted likelihood's own constant),
# offset, (1 | g), (z | g), s() is left to the tape's best of starts, REML, se(), fixed sigma, a prior, a
# bound. Args: arm seeds. Seed 57000 + 17 seed + i.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; seeds <- eval(parse(text = args[2]))
rev_libs(arm)
nug <- if (arm == "lane") 1e-12 else 1e-6
# Vfun(par) gives the covariance without sigma; npar parameters; X fixed effects; w weights; reml
dense <- function(y, X, Vfun, npar, w = NULL, reml = FALSE, se2 = NULL, sig = TRUE, starts) {
  n <- length(y)
  ww <- if (is.null(w)) rep(1, n) else w
  f <- function(p) {
    ls <- if (sig) p[1] else -Inf
    q <- if (sig) p[-1] else p
    V <- Vfun(q)
    diag(V) <- diag(V) + exp(2 * ls) / ww + (if (is.null(se2)) 0 else se2)
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e300)
    wy <- backsolve(R, y, transpose = TRUE); wX <- backsolve(R, X, transpose = TRUE)
    XtX <- crossprod(wX)
    b <- tryCatch(solve(XtX, crossprod(wX, wy)), error = function(e) NULL)
    if (is.null(b)) return(1e300)
    e <- wy - wX %*% b
    v <- 0.5 * sum(e^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
    # a weighted likelihood is not this density: w log N(y | eta, s^2) = log N(y | eta, s^2 / w) + const(s, w)
    if (!is.null(w)) v <- v + 0.5 * sum(ww - 1) * (2 * ls + log(2 * pi)) + 0.5 * sum(log(ww))
    if (reml) v <- v + 0.5 * as.numeric(determinant(XtX)$modulus) - 0.5 * ncol(X) * log(2 * pi)
    if (is.finite(v)) v else 1e300
  }
  best <- NULL
  for (st in starts) {
    o <- tryCatch(nlminb(st, f, control = list(iter.max = 500, eval.max = 900, rel.tol = 1e-12)), error = function(e) NULL)
    if (is.null(o)) next
    o2 <- tryCatch(optim(o$par, f, method = "Nelder-Mead", control = list(maxit = 600, reltol = 1e-13)), error = function(e) NULL)
    if (!is.null(o2) && o2$value < o$objective) { o$par <- o2$par; o$objective <- o2$value }
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  best
}
ke <- function(x, ls, lr) { K <- exp(2 * ls) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * lr))); diag(K) <- diag(K) * (1 + nug); K }
starts_for <- function(np, sg0, sy, rhos = c(0.15, 0.6, 2.5)) {
  out <- list()
  for (r in rhos) for (ns in c(0.3, 0.01)) out[[length(out) + 1L]] <- c(log(sy * ns), rep(c(log(sy), log(r)), length.out = np))
  out
}
res <- list()
for (seed in seeds) for (sg in c(0.1, 1e-2, 1e-3)) {
  set.seed(57000 + 17 * seed + round(-log10(sg)))
  n <- 70L
  x <- sort(round(runif(n, 0, 10), 3)); z <- round(runif(n, 0, 10), 3)
  g <- factor(sample(letters[1:6], n, TRUE)); f2 <- factor(sample(c("a", "b"), n, TRUE))
  ws <- runif(n, 0.5, 2); off <- 0.3 * x
  e <- sg * rnorm(n)
  fx <- 2 * sin(3 * x) + x; fz <- sin(z) * 0.7
  X1 <- matrix(1, n, 1)
  cases <- list(
    two = list(d = data.frame(x, z, y = fx + fz + e), fml = y ~ gp(x) + gp(z), X = X1, np = 4,
               V = function(q) ke(x, q[1], q[2]) + ke(z, q[3], q[4])),
    byf = list(d = data.frame(x, f2, y = fx * (f2 == "a") + sin(x) * (f2 == "b") + e), fml = y ~ gp(x, by = f2), X = cbind(1, as.numeric(f2 == "b")), np = 4,
               V = function(q) { a <- f2 == "a"; K <- matrix(0, n, n); K[a, a] <- ke(x[a], q[1], q[2]); K[!a, !a] <- ke(x[!a], q[3], q[4]); K }),
    byw = list(d = data.frame(x, ws, y = fx * ws + e), fml = y ~ gp(x, by = ws), X = X1, np = 2,
               V = function(q) ke(x, q[1], q[2]) * outer(ws, ws)),
    wts = list(d = data.frame(x, ws, y = fx + e), fml = y | weights(ws) ~ gp(x), X = X1, np = 2, w = ws,
               V = function(q) ke(x, q[1], q[2])),
    off = list(d = data.frame(x, off, y = fx + off + e), fml = y ~ gp(x) + offset(off), X = X1, np = 2, yadj = off,
               V = function(q) ke(x, q[1], q[2])),
    re = list(d = data.frame(x, g, y = fx + rnorm(6, sd = .5)[as.integer(g)] + e), fml = y ~ gp(x) + (1 | g), X = X1, np = 3,
              V = function(q) ke(x, q[1], q[2]) + exp(2 * q[3]) * outer(g, g, "==")),
    reml = list(d = data.frame(x, z, y = fx + 0.2 * z + e), fml = y ~ z + gp(x), X = cbind(1, z), np = 2, reml = TRUE,
                V = function(q) ke(x, q[1], q[2])),
    se = list(d = data.frame(x, s = sg + 1e-4, y = fx + e), fml = y | se(s) ~ gp(x), X = X1, np = 2, se2 = rep((sg + 1e-4)^2, n), sig = FALSE,
              V = function(q) ke(x, q[1], q[2])),
    sesig = list(d = data.frame(x, s = 0.05, y = fx + e + 0.05 * rnorm(n)), fml = y | se(s, sigma = TRUE) ~ gp(x), X = X1, np = 2, se2 = rep(0.05^2, n),
                 V = function(q) ke(x, q[1], q[2]))
  )
  for (cn in names(cases)) {
    cc <- cases[[cn]]
    t0 <- proc.time()[[3]]
    r <- rev_fit(cc$fml, cc$d, REML = isTRUE(cc$reml))
    el <- proc.time()[[3]] - t0
    if (inherits(r$fit, "error")) { cat(sprintf("RST %s %s s=%g seed %d ERROR %s\n", arm, cn, sg, seed, substr(conditionMessage(r$fit), 1, 110))); next }
    ll <- -r$fit$opt$objective
    y0 <- cc$d$y - (if (is.null(cc$yadj)) 0 else cc$yadj)
    sy <- sd(y0)
    sig <- is.null(cc$sig) || cc$sig
    st <- starts_for(cc$np, sg, sy)
    if (cn == "re") st <- lapply(st, function(s) c(s[1:3], log(sy) - 1))
    if (!sig) st <- lapply(st, function(s) s[-1])
    pt <- c(if (sig) log(sigma(r$fit)[1]), r$fit$estimates$theta)
    if (identical(cn, "re")) { ix <- vapply(r$fit$frame$re_blocks, function(b) b$covstruct, ""); th <- r$fit$estimates$theta; gi <- which(ix == "gp"); ui <- which(ix == "us"); pt <- c(if (sig) log(sigma(r$fit)[1]), th[r$fit$frame$re_blocks[[gi]]$theta_idx], th[r$fit$frame$re_blocks[[ui]]$theta_idx]) }
    if (length(pt) == cc$np + sig) st <- c(list(pt), st)
    m <- dense(y0, cc$X, cc$V, cc$np, w = cc$w, reml = isTRUE(cc$reml), se2 = cc$se2, sig = sig, starts = st)
    gap <- ll + m$objective
    res[[cn]] <- rbind(res[[cn]], c(gap, r$fit$opt$convergence, length(r$warn)))
    cat(sprintf("RST %s %-5s s=%g seed %d ll %.4f cf %.4f gap %.2e conv %d w %d %.2fs\n", arm, cn, sg, seed, ll, -m$objective, gap, r$fit$opt$convergence, length(r$warn), el))
    flush.console()
  }
}
for (k in names(res)) { a <- res[[k]]; cat(sprintf("RSTSUM %s %-5s rows %d short>1e-3 %d silent %d worst %.3g over+ %.3g\n", arm, k, nrow(a), sum(a[, 1] < -1e-3), sum(a[, 1] < -1e-3 & a[, 2] == 0 & a[, 3] == 0), min(a[, 1]), max(a[, 1]))) }
# gate edges: a fixed sigma, a prior, a bound
set.seed(58001); x <- sort(runif(60, 0, 10)); d <- data.frame(x, y = 2 * sin(3 * x) + x + 0.1 * rnorm(60))
ed <- function(lab, e) { r <- tryCatch(e, error = function(err) err); if (inherits(r$fit %||% r, "error")) { cat("RED", lab, "ERROR", substr(conditionMessage(r$fit %||% r), 1, 150), "\n"); return(invisible()) }; f <- r$fit; cat(sprintf("RED %s %-22s ll %.4f conv %d w %d sigma %s theta %s\n", arm, lab, -f$opt$objective, f$opt$convergence, length(r$warn), paste(signif(tryCatch(sigma(f)[1], error = function(e) NA), 4)), paste(round(f$estimates$theta, 3), collapse = ","))) }
ed("plain", rev_fit(y ~ gp(x), d))
ed("sigma constant prior", rev_fit(y ~ gp(x), d, prior = set_prior("constant(0.1)", class = "sigma")))
ed("MAP sdgp prior", rev_fit(y ~ gp(x), d, prior = set_prior("normal(0, 1)", class = "sdgp")))
ed("MAP lscale prior", rev_fit(y ~ gp(x), d, prior = set_prior("inv_gamma(3, 1)", class = "lscale")))
ed("bound sigma ub .05", rev_fit(y ~ gp(x), d, prior = set_prior("", class = "sigma", ub = 0.05)))
ed("bound sigma lb 1", rev_fit(y ~ gp(x), d, prior = set_prior("", class = "sigma", lb = 1)))
ed("start beta only", rev_fit(y ~ gp(x), d, start = list(beta = 5)))
ed("start betad only", rev_fit(y ~ gp(x), d, start = list(betad = log(0.5))))
ed("profile=TRUE", rev_fit(y ~ gp(x), d, control = frmtmb_control(profile = TRUE)))
ed("constant X none", rev_fit(y ~ 0 + gp(x), d))
cat("DONE\n")
