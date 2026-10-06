# Punch round 1, m4: the cost of each step of one 300-position kriging
# draw, 300 times, in base R at the sizes of dev/gpby-p1-krigdraw.R.
set.seed(1)
x <- seq(6.05, 9, length.out = 300)
p <- sort(stats::runif(60, 0, 6))
l <- 1.2
D2 <- outer(x, x, "-")^2
Ks <- exp(-outer(x, p, "-")^2 / (2 * l^2))
K <- exp(-outer(p, p, "-")^2 / (2 * l^2)) + diag(1e-6, 60)
Xw <- t(solve(K, t(Ks)))
E <- exp(-D2 / (2 * l^2))
Tm <- tcrossprod(Xw, Ks)
S <- E - Tm
diag(S) <- pmax(1 + 1e-6 - rowSums(Xw * Ks), 0)
R <- chol(S)
tm <- function(lab, f) {
  t <- system.time(for (i in 1:300) f())[["elapsed"]]
  cat(sprintf("%-30s %.3f s for 300\n", lab, t))
}
tm("outer(x, x)^2", function() outer(x, x, "-")^2)
tm("D2 / (2 l^2)", function() D2 / (2 * l^2))
tm("exp(-Q)", function() exp(-D2))
tm("tcrossprod(Xw, Ks)", function() tcrossprod(Xw, Ks))
tm("E - T", function() E - Tm)
tm("chol(S)", function() chol(S))
tm("crossprod(R, z)", function() crossprod(R, stats::rnorm(300)))
tm("solve(K, t(Ks))", function() solve(K, t(Ks)))
tm("Ks kernel 300 x 60", function() exp(-outer(x, p, "-")^2 / (2 * l^2)))
print(extSoftVersion()["BLAS"])
print(La_library())
