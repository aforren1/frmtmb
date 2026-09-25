# multinomial at newdata: trials from newdata, mean counts vs fitted().
source("dev/simnewdata-review/rv-prelude.R")
set.seed(3)
dm <- data.frame(x = rnorm(200), N = sample(5:15, 200, TRUE))
dm$Y <- t(sapply(seq_len(200), function(i) {
  stats::rmultinom(1, dm$N[i], exp(c(0, 0.5 * dm$x[i], -0.3)))
}))
colnames(dm$Y) <- c("a", "b", "c")
fm <- frm(bf(Y | trials(N) ~ x), family = multinomial(K = 3), data = dm)
a <- simulate(fm, nsim = 3, seed = 4)
b <- simulate(fm, nsim = 3, seed = 4, newdata = dm)
cat("identity at fitted rows:", isTRUE(all.equal(a, b, check.attributes = FALSE)), "\n")
ndm <- data.frame(x = c(-1, 1), N = c(2L, 100L))
R <- 4000
sm <- simulate(fm, nsim = R, seed = 1, newdata = ndm)
tots <- vapply(sm, function(v) rowSums(as.matrix(v)), numeric(2))
cat("row totals: row1", paste(range(tots[1, ]), collapse = "-"),
    " row2", paste(range(tots[2, ]), collapse = "-"), "\n")
mc <- Reduce(`+`, lapply(sm, as.matrix)) / R
fvm <- fitted(fm, newdata = ndm)
ex <- fvm[, "Estimate", ]
cat("mean counts:\n"); print(round(mc, 3))
cat("fitted(newdata):\n"); print(round(ex, 3))
