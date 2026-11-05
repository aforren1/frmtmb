# The osa residuals of binary and other discrete families: latent fits against a
# random-effect control on the same build, and the mv cells the misc run missed.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
sumr <- function(r) sprintf("mean %6.3f sd %5.3f min %6.2f max %6.2f", mean(r), sd(r), min(r), max(r))
for (seed in 1:3) {
  set.seed(seed)
  d <- expand.grid(t = 1:10, g = factor(1:15))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(15, arima.sim(list(ar = .5), 10, sd = .6)))
  d$id <- factor(seq_len(n))
  d$ft <- factor(d$t)
  d$bin <- rbinom(n, 1, plogis(.2 + .4 * d$x + e))
  d$bi <- rbinom(n, 5, plogis(.2 + .4 * d$x + e))
  d$tr <- 5L
  d$cnt <- rpois(n, exp(.5 + .4 * d$x + e))
  d$nb <- rnbinom(n, size = 3, mu = exp(.5 + .4 * d$x + e))
  d$ord <- factor(findInterval(.3 * d$x + e + rlogis(n), c(-1, 0, 1)) + 1L, ordered = TRUE)
  specs <- list(
    bern_latent = list(bin ~ x + ar(t, g, cov = TRUE), bernoulli()),
    bern_re = list(bin ~ x + (1 | g), bernoulli()),
    bern_ar1re = list(bin ~ x + ar1(ft + 0 | g), bernoulli()),
    bern_plain = list(bin ~ x, bernoulli()),
    binom_latent = list(bi | trials(tr) ~ x + ar(t, g, cov = TRUE), binomial()),
    binom_re = list(bi | trials(tr) ~ x + (1 | g), binomial()),
    pois_latent = list(cnt ~ x + ar(t, g, cov = TRUE), poisson()),
    nb_latent = list(nb ~ x + ar(t, g, cov = TRUE), negbinomial()),
    ord_latent = list(ord ~ x + ar(t, g, cov = TRUE), cumulative()),
    ord_re = list(ord ~ x + (1 | g), cumulative()))
  for (nm in names(specs)) {
    f <- suppressWarnings(frm(specs[[nm]][[1]], family = specs[[nm]][[2]], data = d))
    r <- msg_of(residuals(f, type = "osa"))
    cat(sprintf("seed %d %-14s %s %s\n", seed, nm, if (is.null(r$err)) sumr(as.matrix(r$val)[, 1]) else paste("ERR", substr(r$err, 1, 110)), if (length(r$warn)) paste("WARN", substr(r$warn[1], 1, 90)) else ""))
  }
}
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:12))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(12, arima.sim(list(ar = .5), 5, sd = .5)))
d$y <- rpois(n, exp(.5 + .4 * d$x + e))
d$pos <- rgamma(n, 4, 4 / exp(.5 + .4 * d$x + e))
d$gs <- .4 * d$x + e + rnorm(n, 0, .5)
d$sub <- d$t != 3
m <- msg_of(frm(bf(y | subset(sub) ~ x + ar(t, g, cov = TRUE)) + bf(gs ~ x), data = d, family = poisson()))
cat("mv subset:", if (is.null(m$err)) as.numeric(logLik(m$val)) else m$err, "\n")
m0 <- msg_of(frm(y | subset(sub) ~ x + ar(t, g, cov = TRUE), data = d, family = poisson()))
m1 <- msg_of(frm(gs ~ x, data = d))
cat("independent sum:", as.numeric(logLik(m0$val)) + as.numeric(logLik(m1$val)), "\n")
m <- msg_of(frm(bf(y ~ x + ar(t, g, cov = TRUE)) + bf(pos ~ x + cosy(t, g), family = Gamma(link = "log")), data = d, family = poisson()))
cat("mv two latents:", if (is.null(m$err)) as.numeric(logLik(m$val)) else m$err, "\n")
if (is.null(m$err)) print(variables(m$val))
