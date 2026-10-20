# Review of lane illse2: step dependence of the route's Hessian on the
# nonlinear body of dev/illse2-rev-asym.R (A = 1, seed 1), and the same
# polynomial written as a linear formula.
#   ILLSE2_ARM=lane Rscript dev/illse2-rev-asymdbg2.R <nl|lin> <sigma>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
kind <- a[1]
sg <- as.numeric(a[2])
A0 <- if (kind == "nl") 1 else 0
deg <- 6L
bn <- paste0("b", 0:deg)
set.seed(1)
n <- 2000
d <- data.frame(x = stats::runif(n, 1, 2), z = stats::runif(n, 0, 10))
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
d$y <- sin(3 * d$x) + A0 * exp(-0.5 * d$z) + stats::rnorm(n, 0, sg)
pf <- paste0("x", 1:deg, collapse = " + ")
l0 <- stats::lm(stats::as.formula(paste("y ~", pf)), data = d)
X <- stats::model.matrix(stats::as.formula(paste("~", pf)), d)
Xm <- mM(X); zm <- m(d$z); ym <- m(d$y)
nb <- deg + 1L
if (kind == "nl") {
  body <- paste(c("b0", paste0("b", 1:deg, " * x", 1:deg),
                  "A * exp(-k * z)"), collapse = " + ")
  fo <- bf(stats::as.formula(paste("y ~", body)),
           stats::as.formula(paste(paste(c(bn, "A", "k"), collapse = " + "),
                                   "~ 1")), nl = TRUE)
  st <- list(beta = c(unname(stats::coef(l0)), A0, 0.5))
  f <- suppressWarnings(frm(fo, data = d, start = st))
  nll <- function(th) {
    mu <- as.vector(Xm %*% col1(th[seq_len(nb)])) +
      th[nb + 1L] * exp(-th[nb + 2L] * zm)
    s2 <- exp(2 * th[nb + 3L])
    (n * log(s2) + sum((ym - mu)^2) / s2) / 2
  }
} else {
  f <- suppressWarnings(frm(stats::as.formula(paste("y ~", pf)), data = d))
  nll <- gauss_nll(X, d$y)
}
p <- f$opt$par
He <- mp_hessian(nll, p)
Ve <- mp_inv(He)
dg <- sqrt(num(diag(Ve)))
V <- suppressWarnings(vcov(f, full = TRUE))
cat(sprintf("%s sigma %g | root %s | vcov cor err %.2e\n", kind, sg,
            !is.null(ns$sdr_of(f)$cov_fixed_root),
            max(abs(num(mM(V) - Ve)) / outer(dg, dg))))
M <- ns$reparam_basis(f$frame, names(p), est = f$estimates)
h <- ns$reparam_hessian_refined(f$obj, p, M)
Mh <- h$M
HMe <- num(crossprod(mM(Mh), He %*% mM(Mh)))
D <- sqrt(diag(HMe))
for (s in c(1e-4, 1e-3, 1e-2, 3e-2, 1e-1)) for (ra in c(FALSE, TRUE)) {
  hh <- ns$reparam_hessian(f$obj, p, Mh, step = s, rich_all = ra)
  Vr <- Mh %*% solve(hh$H) %*% t(Mh)
  cat(sprintf("  step %-6g rich_all %-5s | H err (unit diag) %.2e | asym %.2e | cov cor err %.2e\n",
              s, ra, max(abs(hh$H - HMe) / outer(D, D)),
              max(hh$E) / max(abs(hh$H)),
              max(abs(num(mM(Vr) - Ve)) / outer(dg, dg))))
}
