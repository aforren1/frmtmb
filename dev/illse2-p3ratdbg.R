# Lane illse2, item 3: the rational body of dev/illse2-p3nl.R, seed 3,
# where the lane's route takes a code-1 stop and its covariance is 0.32
# from the exact one on the scale of the correlation. Is the stop the
# maximum? Where does the error sit?
#   ILLSE2_ARM=lane Rscript dev/illse2-p3ratdbg.R [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
seed <- if (length(a)) as.integer(a[1]) else 3L
set.seed(seed)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- sin(3 * d$x) / (1 + 0.3 * d$x) + stats::rnorm(200, 0, 0.3)
f <- bf(y ~ (b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 +
               b5 * x^5 + b6 * x^6) / (1 + dd * x),
        b0 + b1 + b2 + b3 + b4 + b5 + b6 + dd ~ 1, nl = TRUE)
X <- outer(d$x, 0:6, `^`)
Xm <- mM(X); xm <- m(d$x); ym <- m(d$y)
nll <- function(th) {
  mu <- (Xm %*% col1(th[1:7])) / (1 + th[8] * xm)
  r <- ym - mu
  s2 <- exp(2 * th[9])
  (length(d$y) * log(s2) + sum(r^2) / s2) / 2
}
r <- cap(frm(f, data = d))
fit <- r$v
p <- fit$opt$par
for (w in r$w) cat("WARN:", substr(w, 1, 300), "\n")
cat("code", fit$opt$convergence, "| se_explained",
    format(fit$cache$se_explained), "| lost", lost_of(fit), "\n")
He <- mp_hessian(nll, p)
Ve <- mp_inv(He)
ge <- mp_grad(nll, p)
cat("exact Newton gain at the stop:",
    num(sum(ge * (Ve %*% col1(ge)))) / 2, "\n")
D <- sqrt(num(diag(He)))
ev <- eigen(num(He) / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
cat("exact unit-diagonal eigenvalues (double):", signif(ev, 3), "\n")
V <- suppressWarnings(vcov(fit, full = TRUE))
dg <- sqrt(abs(num(diag(Ve))))
cat("SE / exact:", signif(sqrt(diag(V)) / dg, 4), "\n")
# a restart: is there a better point?
o <- stats::optim(p, function(t) num(nll(m(t))), method = "BFGS",
                  control = list(reltol = 1e-14, maxit = 2000))
cat("BFGS on the exact nll from the stop gains", num(nll(m(p))) - o$value,
    "\n")
# the route's Hessian under both difference schemes, against the exact
M <- ns$reparam_basis(fit$frame, names(p), est = fit$estimates)
cat("M covers", length(attr(M, "pos")), "positions\n")
Mx <- mM(M)
Hmx <- t(Mx) %*% He %*% Mx
for (ra in if (nzchar(Sys.getenv("I2_SCHEMES"))) c(FALSE, TRUE)) {
  h <- ns$reparam_hessian(fit$obj, p, M, rich_all = ra)
  Dm <- sqrt(diag(h$H))
  em <- eigen(h$H / outer(Dm, Dm), symmetric = TRUE, only.values = TRUE)$values
  rel <- max(abs(num(mM(h$H) - Hmx))) / max(abs(num(Hmx)))
  Vr <- M %*% solve(h$H) %*% t(M)
  cat(sprintf("rich_all %s | M-kappa %.2e | H_M entry error %.2e | SE / exact %s\n",
              ra, max(em) / min(em), rel,
              paste(signif(sqrt(diag(Vr)) / dg, 4), collapse = " ")))
}
# the natural inverses: of the exact AD Hessian, and the one vcov() shows
Vad <- solve(fit$obj$he(p))
Vsh <- suppressWarnings(vcov(fit, full = TRUE))
ce <- function(V) max(abs(num(mM(V) - Ve)) / outer(dg, dg))
cat(sprintf("cor err: solve(AD Hessian) %.2e | vcov(full) %.2e\n",
            ce(Vad), ce(Vsh)))
