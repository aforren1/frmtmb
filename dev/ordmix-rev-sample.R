# Reviewer of lane ordmix: frmtmb.sample on an ordinal mixture the lane's
# draws test does not cover: theta1 ~ z and cs() on component 2 only.
# posterior_epred() and log_lik() at stored draws against a hand
# computation from the draws' brms-named columns; posterior_predict()
# frequencies against posterior_epred(). Data seed 20261098, sampler
# seed 5, one chain of 400.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", find.package("frmtmb"), " sample",
    find.package("frmtmb.sample"), "\n")
set.seed(20261098)
n <- 250
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- rbinom(n, 1, plogis(0.4 * d$z))
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
st <- "student_t(3, 0, 2.5)"
pr <- set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
  set_prior("normal(0, 2)", class = "b", dpar = "mu2") +
  set_prior("normal(0, 2)", class = "b", dpar = "theta1") +
  set_prior(st, class = "Intercept", dpar = "mu1") +
  set_prior(st, class = "Intercept", dpar = "mu2")
fit <- frm(bf(y ~ x, theta1 ~ z, mu2 ~ x + cs(z)),
           family = mixture(cumulative(), sratio()), data = d, prior = pr)
print(fixef(fit)[, 1:2])
ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 5)
M <- as_draws_matrix(ds)
cat("draw columns:", paste(colnames(M), collapse = " "), "\n")
E <- posterior_epred(ds)
L <- log_lik(ds)
cat("epred dims", paste(dim(E), collapse = "x"), " log_lik dims",
    paste(dim(L), collapse = "x"), "\n")
hand <- function(s) {
  m <- setNames(as.numeric(M[s, ]), colnames(M))
  t1 <- m[paste0("b_mu1_Intercept[", 1:3, "]")]
  t2 <- m[paste0("b_mu2_Intercept[", 1:3, "]")]
  cs <- m[paste0("bcs_mu2_z[", 1:3, "]")]
  e1 <- m[["b_mu1_x"]] * d$x
  e2 <- m[["b_mu2_x"]] * d$x
  th <- plogis(m[["b_theta1_Intercept"]] + m[["b_theta1_z"]] * d$z)
  p1 <- t(sapply(seq_len(n), function(i) diff(c(0, plogis(t1 - e1[i]), 1))))
  p2 <- t(sapply(seq_len(n), function(i) {
    h <- plogis(t2 - cs * d$z[i] - e2[i])
    c(h[1], (1 - h[1]) * h[2], (1 - h[1]) * (1 - h[2]) * h[3], prod(1 - h))
  }))
  th * p1 + (1 - th) * p2
}
for (s in c(1, 100, nrow(M))) {
  P <- hand(s)
  ll <- log(P[cbind(seq_len(n), d$y)])
  cat(sprintf("draw %d: max|epred - hand| %.3g, max|log_lik - hand| %.3g\n",
              s, max(abs(E[s, , ] - P)), max(abs(L[s, ] - ll))))
}
Y <- posterior_predict(ds)
pb <- colMeans(apply(E, c(2, 3), mean))
obs <- sapply(1:4, function(k) mean(Y == k))
cat("posterior_predict category shares:", sprintf("%.4f", obs), "\n")
cat("posterior_epred   category shares:", sprintf("%.4f", pb), "\n")
# a binomial standard error on draws x rows cells, an overstatement of
# independence that makes the check conservative in the other direction
se <- sqrt(pb * (1 - pb) / length(Y))
cat("z:", sprintf("%.2f", (obs - pb) / se), "\n")
