# Attack 7: refusals, against brms run on each call, and the guards' absent
# cases. A brms "build" here needs make_stancode() AND a compile through
# rstan::stanc() (the parse and the typecheck of the program); the full C++
# compile of a sampled set is in the second half.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(brms); library(rstan)})
set.seed(20261102)
d <- expand.grid(t = 1:5, g = factor(1:10))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- rpois(n, exp(0.5 + 0.3 * d$x))
d$bin <- rbinom(n, 1, 0.5)
d$tr <- 5L
d$bi <- rbinom(n, 5, 0.5)
d$pos <- rgamma(n, 3, 3)
d$un <- rbeta(n, 2, 2)
d$ord <- factor(sample(1:3, n, TRUE), ordered = TRUE)
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$gs <- rnorm(n)
d$id <- seq_len(n)
d$id2 <- factor(seq_len(n))
M <- diag(n)
data2 <- list(V = M, W = {
  W <- matrix(0, n, n)
  for (i in 1:(n - 1)) W[i, i + 1] <- W[i + 1, i] <- 1
  W / rowSums(W)
})
br <- function(form_q, fam_q) {
  ns <- asNamespace("brms")
  r <- msg_of(brms::stancode(eval(form_q, list(), ns), d, family = eval(fam_q, list(), ns), data2 = data2))
  if (!is.null(r$err)) return(c("REFUSES", substr(gsub("\\s+", " ", r$err), 1, 170)))
  chk <- msg_of(rstan::stanc(model_code = r$val, allow_undefined = FALSE, verbose = FALSE))
  if (!is.null(chk$err)) return(c("BUILDS,STANC-FAILS", substr(chk$err, 1, 170)))
  c("BUILDS", "")
}
fr <- function(form_q, fam_q) {
  ns <- asNamespace("frmtmb")
  r <- msg_of(frm(eval(form_q, list(), ns), data = d, family = eval(fam_q, list(), ns), data2 = data2, dry_run = "frame"))
  if (!is.null(r$err)) return(c("REFUSES", substr(gsub("\\s+", " ", r$err), 1, 170)))
  c("BUILDS", "")
}
suppressPackageStartupMessages({library(frmtmb.spatial); library(frmtmb.survival); library(frmtmb.eam)})
cases <- list(
  list("ma cov=F poisson", quote(y ~ x + ma(t, g)), quote(poisson())),
  list("arma cov=F poisson", quote(y ~ x + arma(t, g)), quote(poisson())),
  list("ar p2 cov=T poisson", quote(y ~ x + ar(t, g, p = 2, cov = TRUE)), quote(poisson())),
  list("ma q2 cov=T poisson", quote(y ~ x + ma(t, g, q = 2, cov = TRUE)), quote(poisson())),
  list("arma p2 cov=T Gamma", quote(pos ~ x + arma(t, g, p = 2, cov = TRUE)), quote(Gamma(link = "log"))),
  list("fcor poisson", quote(y ~ x + fcor(V)), quote(poisson())),
  list("sar poisson", quote(y ~ x + sar(W)), quote(poisson())),
  list("sar error poisson", quote(y ~ x + sar(W, type = "error")), quote(poisson())),
  list("sar bernoulli", quote(bin ~ x + sar(W)), quote(bernoulli())),
  list("ar categorical", quote(cat ~ x + ar(t, g, cov = TRUE)), quote(categorical())),
  list("ar multinomial", quote(bf(cbind(bi, tr) | trials(tr) ~ x + ar(t, g, cov = TRUE))), quote(multinomial())),
  list("ar dirichlet", quote(bf(cbind(un, 1 - un) ~ x + ar(t, g, cov = TRUE))), quote(dirichlet())),
  list("ar mixture", quote(y ~ x + ar(t, g, cov = TRUE)), quote(mixture(poisson(), poisson()))),
  list("ar mixture gauss", quote(gs ~ x + ar(t, g, cov = TRUE)), quote(mixture(gaussian(), gaussian()))),
  list("ar on shape (dpar)", quote(bf(y ~ x, shape ~ ar(t, g, cov = TRUE))), quote(negbinomial())),
  # the nl bodies and a nonlinear dpar beside a linear mu: must NOT fire
  list("ar on mu, nlf(shape)", quote(bf(y ~ x + ar(t, g, cov = TRUE), nlf(shape ~ exp(a)), a ~ 1)), quote(negbinomial())),
  list("ar on mu, shape linear", quote(bf(y ~ x + ar(t, g, cov = TRUE), shape ~ x)), quote(negbinomial())),
  list("ar on mu, zi nl", quote(bf(y ~ x + ar(t, g, cov = TRUE), nlf(zi ~ inv_logit(a)), a ~ 1)), quote(zero_inflated_poisson())),
  list("nl mu + ar", quote(bf(y ~ a * exp(b * x) + ar(t, g, cov = TRUE), a + b ~ 1, nl = TRUE)), quote(poisson())),
  list("nl mu + autocor", quote(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE, autocor = ~ ar(t, g, cov = TRUE))), quote(poisson())),
  # guards that must not fire on a correct model
  list("ar cov=F p2 poisson", quote(y ~ x + ar(t, g, p = 2)), quote(poisson())),
  list("ar cov=T p1 poisson", quote(y ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("ma cov=T q1 poisson", quote(y ~ x + ma(t, g, cov = TRUE)), quote(poisson())),
  list("arma 1,1 cov=T bernoulli", quote(bin ~ x + arma(t, g, cov = TRUE)), quote(bernoulli())),
  list("ar cov=T p2 gaussian", quote(gs ~ x + ar(t, g, p = 2, cov = TRUE)), quote(gaussian())),
  list("ma cov=F gaussian", quote(gs ~ x + ma(t, g)), quote(gaussian())),
  list("ma cov=F student", quote(gs ~ x + ma(t, g)), quote(student())),
  list("fcor gaussian", quote(gs ~ x + fcor(V)), quote(gaussian())),
  # other families and terms the grid did not name
  list("von_mises ar", quote(I(gs/10) ~ x + ar(t, g, cov = TRUE)), quote(von_mises())),
  list("cox ar", quote(pos ~ x + ar(t, g, cov = TRUE)), quote(cox())),
  list("asym_laplace ar", quote(gs ~ x + ar(t, g, cov = TRUE)), quote(asym_laplace())),
  list("shifted_lognormal ar", quote(pos ~ x + ar(t, g, cov = TRUE)), quote(shifted_lognormal())),
  list("hurdle_gamma ar", quote(pos ~ x + ar(t, g)), quote(hurdle_gamma())),
  list("hurdle_lognormal ar", quote(pos ~ x + ar(t, g)), quote(hurdle_lognormal())),
  list("zero_inflated_beta ar", quote(un ~ x + ar(t, g)), quote(zero_inflated_beta())),
  list("zero_inflated_negbinomial ar", quote(y ~ x + ar(t, g)), quote(zero_inflated_negbinomial())),
  list("zero_inflated_binomial ar", quote(bi | trials(tr) ~ x + ar(t, g)), quote(zero_inflated_binomial())),
  list("gen_extreme_value ar", quote(gs ~ x + ar(t, g, cov = TRUE)), quote(gen_extreme_value())),
  list("wiener ar", quote(pos | dec(bin) ~ x + ar(t, g, cov = TRUE)), quote(wiener())),
  list("xbeta ar", quote(un ~ x + ar(t, g, cov = TRUE)), quote(xbeta())),
  # time and group shape conditions
  list("all single-row groups", quote(y ~ x + ar(id, id2, cov = TRUE)), quote(poisson())),
  list("one row per group cosy", quote(y ~ x + cosy(gr = id2)), quote(poisson()))
)
res <- list()
for (cc in cases) {
  b <- br(cc[[2]], cc[[3]])
  f <- fr(cc[[2]], cc[[3]])
  same <- b[1] == f[1] || (b[1] == "BUILDS,STANC-FAILS" && f[1] == "REFUSES")
  cat(sprintf("%-30s brms %-18s frm %-8s %s\n", cc[[1]], b[1], f[1], if (same) "" else "<<< DIFF"))
  if (nzchar(b[2])) cat("    brms:", b[2], "\n")
  if (nzchar(f[2])) cat("    frm: ", f[2], "\n")
}
