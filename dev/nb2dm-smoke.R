# First look at the two families on the lane build.
#   Rscript dev/nb2dm-smoke.R > dev/nb2dm-out/smoke.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "\n")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(15, 20))
d$y <- rnbinom(n, mu = exp(1 + 0.4 * d$x), size = 2.5)
f1 <- frm(y ~ x, family = negbinomial(), data = d)
f2 <- frm(y ~ x, family = brmsfamily("negbinomial2"), data = d)
print(f2)
cat("logLik nb, nb2:", format(logLik(f1), digits = 15),
    format(logLik(f2), digits = 15), "\n")
print(fixef(f1)); print(fixef(f2))
cat("sigma(f2):", format(sigma(f2), digits = 15), " 1 / shape(f1):",
    format(1 / exp(f1$estimates$betad[[1]]), digits = 15), "\n")
str(f2$estimates)
print(head(fitted(f2), 3)); print(head(fitted(f1), 3))
print(head(predict(f2), 3))
print(head(residuals(f2, type = "pearson"), 3))
print(head(residuals(f1, type = "pearson"), 3))
s <- simulate(f2, nsim = 2, seed = 1); print(dim(s))
print(summary(f2))
f3 <- frm(bf(y ~ x, sigma ~ z), family = brmsfamily("negbinomial2"),
          data = d)
f3b <- frm(bf(y ~ x, shape ~ z), family = negbinomial(), data = d)
cat("logLik sigma~z vs shape~z:", format(logLik(f3), digits = 15),
    format(logLik(f3b), digits = 15), "\n")
print(fixef(f3)); print(fixef(f3b))
d$e <- exp(runif(n))
f4 <- frm(y | rate(e) ~ x, family = brmsfamily("negbinomial2"), data = d)
f4b <- frm(y | rate(e) ~ x, family = negbinomial(), data = d)
cat("logLik rate:", format(logLik(f4), digits = 15),
    format(logLik(f4b), digits = 15), "\n")
f5 <- frm(y ~ x + (1 | g), family = brmsfamily("negbinomial2"), data = d)
print(f5)
print(family(f2)); print(family(f2)$link_sigma)
print(brmsfamily("negbinomial2", link_sigma = "softplus")$link_sigma)
print(try(brmsfamily("negbinomial2", link = "logit")))

## dirichlet multinomial
set.seed(12)
K <- 3
tr <- sample(10:30, n, TRUE)
eta <- cbind(0, 0.3 + 0.5 * d$x, -0.2 - 0.4 * d$x)
P <- exp(eta) / rowSums(exp(eta))
phi <- 8
Y <- t(sapply(seq_len(n), function(i) {
  p <- rgamma(K, P[i, ] * phi); p <- p / sum(p)
  rmultinom(1, tr[i], p)
}))
colnames(Y) <- c("a", "b", "c")
d$Y <- Y; d$tr <- tr
g1 <- frm(Y | trials(tr) ~ x, family = dirichlet_multinomial(), data = d)
print(g1)
g0 <- frm(Y | trials(tr) ~ x, family = multinomial(), data = d)
cat("logLik dm, mn:", format(logLik(g1), digits = 15),
    format(logLik(g0), digits = 15), "\n")
print(fixef(g1)); print(fixef(g0))
print(dim(fitted(g1))); print(fitted(g1)[1:2, , ])
print(dim(predict(g1)))
print(dim(simulate(g1, nsim = 1, seed = 2)[[1]]))
print(family(g1)); print(family(g1)$link_phi)
g2 <- frm(bf(Y | trials(tr) ~ x, phi ~ z), family = dirichlet_multinomial(),
          data = d)
print(fixef(g2))
g3 <- frm(Y | trials(tr) ~ x + (1 | g), family = dirichlet_multinomial(),
          data = d)
print(g3)
g4 <- frm(Y | trials(tr) ~ x, family = dirichlet_multinomial(refcat = "b"),
          data = d)
cat("refcat b logLik:", format(logLik(g4), digits = 15), "\n")
print(try(residuals(g1)))
print(try(dharma_residuals(g1)))
print(frm_compat("dirichlet_multinomial"))
