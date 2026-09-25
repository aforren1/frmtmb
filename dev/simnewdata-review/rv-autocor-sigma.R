# Recheck of rv-autocor.R case (h), whose variance read 2.7% high at
# seed 5: two more seeds, 20000 draws each.
source("dev/simnewdata-review/rv-prelude.R")
set.seed(8)
n_g <- 60
d <- expand.grid(time = 1:6, g = factor(seq_len(n_g)))
d$y <- 1 + unlist(lapply(seq_len(n_g), function(i) {
  as.vector(stats::arima.sim(list(ar = 0.7), 6))
}))
d$x <- stats::rnorm(nrow(d))
d$y3 <- d$y + (exp(0.5 * d$x) - 1) * stats::rnorm(nrow(d))
fits <- frm(bf(y3 ~ 1 + ar(time, gr = g, cov = TRUE), sigma ~ x), data = d)
nds <- data.frame(time = c(2, 3, 4), g = factor("A"), x = c(-1, 0, 1.5))
sgr <- as.vector(frm_linpred(fits, newdata = nds, dpar = "sigma",
                             type = "response"))
S <- outer(sgr, sgr) * autocor_matrix(fits)[2:4, 2:4]
for (sd_ in c(6, 7)) {
  ss <- as.matrix(simulate(fits, nsim = 20000, seed = sd_, newdata = nds))
  E <- stats::cov(t(ss))
  z <- (E - S) / sqrt((outer(diag(S), diag(S)) + S^2) / 20000)
  cat("seed", sd_, "cov z:", format(z[upper.tri(z, TRUE)], digits = 2), "\n")
}
