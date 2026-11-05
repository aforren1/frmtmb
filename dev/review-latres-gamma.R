# Attack 3: dispersion families with one latent per row. Does a plain Gamma
# with data from the model give a sane logLik? Seeds 1..N of the same DGP, both
# forms of the term, against a Gamma with no latent, and the shape at the stop.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
one <- function(seed, form, family, resp) {
  set.seed(seed)
  d <- expand.grid(t = 1:6, g = factor(1:25))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(25, arima.sim(list(ar = 0.5), 6, sd = 0.5)))
  d$pos <- rgamma(n, 4, 4 / exp(0.5 + 0.4 * d$x + e))
  d$nb <- rnbinom(n, size = 4, mu = exp(1.5 + 0.4 * d$x + e))
  m <- msg_of(frm(form, data = d, family = family))
  if (!is.null(m$err)) return(data.frame(seed = seed, ll = NA, conv = NA, b1 = NA, sderr = NA, nat = NA, disp = NA, nwarn = NA, err = substr(m$err, 1, 80)))
  f <- m$val
  th <- f$estimates$theta
  v <- variables(f)
  eap <- tryCatch(summary(f)$spec_pars, error = function(e) NULL)
  disp <- if (!is.null(eap) && nrow(eap)) eap[1, 1] else NA
  sm <- tryCatch(summary(f)$cor_pars, error = function(e) NULL)
  data.frame(seed = seed, ll = as.numeric(logLik(f)), conv = f$opt$convergence,
             b1 = f$estimates$beta[2], sderr = exp(th[1]),
             nat = if (!is.null(sm)) sm[1, 1] else NA, disp = disp, nwarn = length(m$warn),
             err = "")
}
run <- function(label, form, family) {
  r <- do.call(rbind, lapply(1:12, one, form = form, family = family))
  cat("\n==", label, "\n")
  print(r, digits = 5, row.names = FALSE)
}
run("Gamma log, ar cov=FALSE (on response)", pos ~ x + ar(t, g), Gamma(link = "log"))
run("Gamma log, ar cov=TRUE (link scale)", pos ~ x + ar(t, g, cov = TRUE), Gamma(link = "log"))
run("Gamma log, cosy", pos ~ x + cosy(t, g), Gamma(link = "log"))
run("negbinomial, ar cov=TRUE", nb ~ x + ar(t, g, cov = TRUE), negbinomial())
run("negbinomial, ar cov=FALSE (log scale)", nb ~ x + ar(t, g), negbinomial())
run("lognormal, ar cov=TRUE", pos ~ x + ar(t, g, cov = TRUE), lognormal())
