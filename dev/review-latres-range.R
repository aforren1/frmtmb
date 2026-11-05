# Attack 3: latres_mu_range() is a list of family names. Which families with a
# post-link AR part (on_response TRUE) are outside it, and what does the tape say
# when a latent pushes mu off their support? Stan rejects there.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(3)
d <- expand.grid(t = 1:5, g = factor(1:10))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$cnt <- rpois(n, exp(1 + .3 * d$x))
d$pos <- rgamma(n, 4, 4 / exp(.3 * d$x))
d$un <- pmin(pmax(rbeta(n, 3, 3), 1e-3), 1 - 1e-3)
d$bi <- rbinom(n, 8, .4)
d$tr <- 8L
fams <- list(
  list("xbeta", un ~ x + ar(t, g), xbeta()),
  list("discrete_weibull", cnt ~ x + ar(t, g), discrete_weibull()),
  list("nbinom1", cnt ~ x + ar(t, g), frmtmb::nbinom1()),
  list("tweedie", pos ~ x + ar(t, g), frmtmb::tweedie()),
  list("hurdle_cumulative", cnt ~ x + ar(t, g), "hurdle_cumulative"),
  list("beta_binomial", bi | trials(tr) ~ x + ar(t, g), beta_binomial()),
  list("zero_inflated_beta_binomial", bi | trials(tr) ~ x + ar(t, g), zero_inflated_beta_binomial()),
  list("shifted_lognormal", pos ~ x + ar(t, g), shifted_lognormal()),
  list("Gamma log (control)", pos ~ x + ar(t, g), Gamma(link = "log")),
  list("Beta (control)", un ~ x + ar(t, g), Beta()),
  list("poisson identity (control)", cnt ~ x + ar(t, g), poisson("identity")))
for (ff in fams) {
  m <- msg_of(frm(ff[[2]], data = d, family = ff[[3]]))
  if (!is.null(m$err)) {
    cat(sprintf("%-30s ERR %s\n", ff[[1]], substr(m$err, 1, 150)))
    next
  }
  f <- m$val
  bk <- fi("frame_latres")(f$frame)[[1]]
  a <- bk$aux
  env <- f$obj$env
  # a latent that pushes row 2's mean far below zero, or above one
  full <- env$last.par.best
  nm <- names(full)
  full[nm == "theta"][1] <- log(1)
  full[nm == "theta"][2] <- 3
  z <- rep(0, sum(nm == "b"))
  z[1] <- 8
  full[nm == "b"] <- z
  v1 <- env$f(full)
  z[1] <- -8
  full[nm == "b"] <- z
  v2 <- env$f(full)
  cat(sprintf("%-30s on_response=%-5s range=%-5s code %d  nll at z1=+8: %10.3f  z1=-8: %10.3f  ll %.3f\n",
              ff[[1]], isTRUE(a$on_response), a$mu_range %||% "none", f$opt$convergence, v1, v2, as.numeric(logLik(f))))
}
