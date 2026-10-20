# Review of lane illse2: why the nonlinear body of dev/illse2-rev-asym.R
# (A = 1, seed 1) keeps a covariance 3e-4 from the exact one.
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
A0 <- 1; deg <- 6L
bn <- paste0("b", 0:deg)
body <- paste(c("b0", paste0("b", 1:deg, " * x", 1:deg),
                "A * exp(-k * z)"), collapse = " + ")
fo <- bf(stats::as.formula(paste("y ~", body)),
         stats::as.formula(paste(paste(c(bn, "A", "k"), collapse = " + "),
                                 "~ 1")), nl = TRUE)
set.seed(1)
n <- 2000
d <- data.frame(x = stats::runif(n, 1, 2), z = stats::runif(n, 0, 10))
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
d$y <- sin(3 * d$x) + A0 * exp(-0.5 * d$z) + stats::rnorm(n, 0, 0.01)
l0 <- stats::lm(stats::as.formula(paste("y ~", paste0("x", 1:deg,
                                  collapse = " + "))), data = d)
st <- list(beta = c(unname(stats::coef(l0)), A0, 0.5))
f <- suppressWarnings(frm(fo, data = d, start = st))
p <- f$opt$par
print(rownames(fixef(f)))
print(names(p))
X <- stats::model.matrix(stats::as.formula(paste("~", paste0("x", 1:deg,
                                  collapse = " + "))), d)
Xm <- mM(X); zm <- m(d$z); ym <- m(d$y)
nb <- deg + 1L
nll <- function(th) {
  mu <- Xm %*% col1(th[seq_len(nb)])
  mu <- as.vector(mu) + th[nb + 1L] * exp(-th[nb + 2L] * zm)
  s2 <- exp(2 * th[nb + 3L])
  (n * log(s2) + sum((ym - mu)^2) / s2) / 2
}
cat("nll mp - obj$fn:", num(nll(m(p))) - f$obj$fn(p), "\n")
g <- mp_grad(nll, p)
cat("mp grad:", format(num(g), digits = 3), "\n")
cat("AD grad:", format(as.numeric(f$obj$gr(p)), digits = 3), "\n")
He <- mp_hessian(nll, p)
Hd <- f$obj$he(p)
cat("AD Hessian vs exact, max rel entry:", max(abs(Hd - num(He)) /
                                                 max(abs(num(He)))), "\n")
Ve <- mp_inv(He)
dg <- sqrt(num(diag(Ve)))
V <- suppressWarnings(vcov(f, full = TRUE))
E <- abs(num(mM(V) - Ve)) / outer(dg, dg)
cat("cor err by row:", format(apply(E, 1, max), digits = 2), "\n")
cat("SE ratio:", format(sqrt(diag(V)) / dg, digits = 8), "\n")
sdr <- ns$sdr_of(f)
cat("root:", !is.null(sdr$cov_fixed_root), " lost:", lost_of(f), "\n")
M <- ns$reparam_basis(f$frame, names(p), est = f$estimates)
cat("M covers:", attr(M, "pos"), "\n")
h <- ns$reparam_hessian_refined(f$obj, p, M)
Mh <- h$M
HMe <- crossprod(mM(Mh), He %*% mM(Mh))
D <- sqrt(num(diag(HMe)))
cat("route H in its coordinates vs exact M'HM (unit diag), max:",
    max(abs(h$H - num(HMe)) / outer(D, D)), "\n")
ev <- eigen(num(HMe) / outer(D, D), only.values = TRUE)$values
cat("kappa of exact M'HM at unit diag:", max(ev) / min(ev), "\n")
cat("asym:", max(h$E) / max(abs(h$H)), "\n")
# the covariance from the exact M'HM, in 256 bits
Vr <- mM(Mh) %*% mp_inv(HMe) %*% t(mM(Mh))
cat("M (M'HM)^-1 M' vs exact, cor:", max(abs(num(Vr - Ve)) / outer(dg, dg)),
    "\n")
Vroute <- Mh %*% solve(h$H) %*% t(Mh)
cat("M solve(route H) M' vs exact, cor:",
    max(abs(num(mM(Vroute) - Ve)) / outer(dg, dg)), "\n")
