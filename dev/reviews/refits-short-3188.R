# Short check of the refits fix round, 3188 LOW: the marginaleffects data
# attribute after refit(). Seed 5, n = 200.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
set.seed(2501)
dm <- data.frame(x = rnorm(60), g2 = rnorm(60))
dm$yd <- sample(c("1", "2", "3", "4"), 60, TRUE)
go("1259 mvbf: character ordinal beside a gaussian",
   frm(frmtmb::bf(yd ~ x) + frmtmb::bf(g2 ~ x), data = dm,
       family = list(frmtmb::cumulative(), gaussian())))
set.seed(5)
n <- 200
x <- rnorm(n)
N <- rep(12L, n)
p <- cbind(a = 1, b = exp(0.3 + 0.4 * x), c = exp(-0.2 * x))
p <- p / rowSums(p)
Y <- t(vapply(seq_len(n), function(i) rmultinom(1, N[i], p[i, ]), numeric(3)))
colnames(Y) <- c("a", "b", "c")
dd <- data.frame(x = x, N = N)
dd$Y <- Y
fm <- frm(Y | trials(N) ~ x, data = dd, family = multinomial())
md0 <- attr(fm, "marginaleffects_modeldata")
cat("attr names:", names(md0), "\n")
set.seed(6)
Y2 <- Y[sample(n), ]
rf <- refit(fm, Y2)
md1 <- attr(rf, "marginaleffects_modeldata")
go("R1 refit attr Y.a..Y.c == new response",
   identical(unname(as.matrix(md1[, c("Y.a", "Y.b", "Y.c")])), unname(Y2 * 1)))
go("R2 refit attr x unchanged", identical(md1$x, md0$x))
go("R3 refit data_frame Y == new response",
   identical(unname(rf$frame$data_frame$Y * 1), unname(Y2 * 1)))
go("R4 original attr unchanged",
   identical(unname(as.matrix(md0[, c("Y.a", "Y.b", "Y.c")])), unname(Y * 1)))
go("R5 refit coef differ from fit", !identical(fixef(rf), fixef(fm)))
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  go("R6 avg_predictions(rf) estimates",
     round(marginaleffects::avg_predictions(rf)$estimate, 6))
  go("R6 colMeans(fitted(rf))",
     round(colMeans(fitted(rf)[, "Estimate", ]), 6))
}
# cbind() spelling: the raw columns are carried beside the matrix
dd2 <- data.frame(x = x, N = N, Ya = Y[, 1], Yb = Y[, 2], Yc = Y[, 3])
fc <- frm(cbind(Ya, Yb, Yc) | trials(N) ~ x, data = dd2, family = multinomial())
mc0 <- attr(fc, "marginaleffects_modeldata")
go("R7 cbind fit attr names", names(mc0))
rc <- refit(fc, Y2)
mc1 <- attr(rc, "marginaleffects_modeldata")
go("R8 cbind refit attr names", names(mc1))
go("R9 cbind refit: attr Ya == new response column 1",
   if (!is.null(mc1) && "Ya" %in% names(mc1)) identical(mc1$Ya, Y2[, 1] * 1) else "no Ya")
go("R10 cbind refit: data_frame response == new", {
  df <- rc$frame$data_frame
  nm <- grep("cbind", names(df), value = TRUE)
  list(nm, if (length(nm)) identical(unname(df[[nm]] * 1), unname(Y2 * 1)))
})
# vcov_cluster(~ v) on a refit of a matrix-response fit, v not a model term
dd$cl <- factor(rep(1:20, each = 10))
fm2 <- frm(Y | trials(N) ~ x, data = dd, family = multinomial())
rf2 <- refit(fm2, Y2)
go("R11 vcov_cluster(fm2, ~ cl)", se <- round(sqrt(diag(vcov_cluster(fm2, ~ cl))), 6))
go("R12 vcov_cluster(refit, ~ cl)", round(sqrt(diag(vcov_cluster(rf2, ~ cl))), 6))
go("R13 vcov_cluster(refit, dd$cl)", round(sqrt(diag(vcov_cluster(rf2, dd$cl))), 6))
# the same for a univariate refit (data_frame not rewritten)
du <- data.frame(x = x, y = rnorm(n), cl = dd$cl)
fu <- frm(y ~ x, data = du)
go("R14 vcov_cluster(refit(univariate), ~ cl)",
   round(sqrt(diag(vcov_cluster(refit(fu, rnorm(n)), ~ cl))), 6))
