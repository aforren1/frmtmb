# Short check of the refits fix round: follow-ups of F2 and the 3188
# refresh. ARM=before reads rellib-r21. Seeds inline.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
G <- 20
mkd <- function() {
  set.seed(1)
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   cl = factor(rep(seq_len(G / 2), each = 12)),
                   x = rnorm(G * 6))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G * 6) + rnorm(G / 2)[dd$cl]
  dd
}
se_of <- function(v) round(sqrt(diag(v)), 7)
ref <- se_of(vcov_cluster(frm(y ~ x, data = mkd()), mkd()$cl))
go("X0 reference (dd$cl)", ref)
mk <- function() { d_local <- mkd(); frm(y ~ x, data = d_local) }
fl <- mk()
cl <- factor(rep(1:2, length.out = 120))
go("X1 data gone (local), global cl of length n: used or refused?",
   se_of(vcov_cluster(fl, ~ cl)))
go("X1 explicit global cl", se_of(vcov_cluster(fl, cl)))
rm(cl)

cat("\n== matrix-response refit and vcov_cluster(~ v) ==\n")
set.seed(5)
n <- 200
x <- rnorm(n)
N <- rep(12L, n)
p <- cbind(a = 1, b = exp(0.3 + 0.4 * x), c = exp(-0.2 * x))
p <- p / rowSums(p)
Y <- t(vapply(seq_len(n), function(i) rmultinom(1, N[i], p[i, ]), numeric(3)))
colnames(Y) <- c("a", "b", "c")
dd <- data.frame(x = x, N = N, cl = factor(rep(1:20, each = 10)))
dd$Y <- Y
fm2 <- frm(Y | trials(N) ~ x, data = dd, family = multinomial())
set.seed(6)
Y2 <- Y[sample(n), ]
rf2 <- refit(fm2, Y2)
go("M1 vcov_cluster(refit, ~ cl)", se_of(vcov_cluster(rf2, ~ cl)))
go("M2 vcov_cluster(refit, dd$cl)", se_of(vcov_cluster(rf2, dd$cl)))

cat("\n== influence() of a univariate refit: which response? ==\n")
set.seed(8)
du <- data.frame(x = rnorm(30), y = rnorm(30))
fu <- frm(y ~ x, data = du)
ynew <- rnorm(30, 2)
ru <- refit(fu, ynew)
iu <- suppressWarnings(influence(ru, force = TRUE))
du2 <- du; du2$y <- ynew
i2 <- suppressWarnings(influence(frm(y ~ x, data = du2), force = TRUE))
i0 <- suppressWarnings(influence(fu, force = TRUE))
go("U1 influence(refit) fixed = influence(fresh fit on ynew)",
   max(abs(as.matrix(iu$fixed) - as.matrix(i2$fixed))))
go("U2 influence(refit) fixed = influence(original fit)",
   max(abs(as.matrix(iu$fixed) - as.matrix(i0$fixed))))
