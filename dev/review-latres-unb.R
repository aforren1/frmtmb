# Short check, item 1: "an out-of-range mean is the only way this block reaches an
# unbounded density" (findings, M4). Held links, no refused form: a dispersion
# that the rows do not ask for, 30 seeds per family.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
run <- function(lab, sim, form, fam, plain) {
  out <- list()
  for (s in 1:30) {
    set.seed(s)
    d <- expand.grid(t = 1:6, g = factor(1:25))
    d <- d[order(d$g, d$t), ]
    n <- nrow(d)
    d$x <- rnorm(n)
    d$e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
    d$tr <- 10L
    d$y <- sim(d)
    w <- character(0)
    f <- withCallingHandlers(tryCatch(frm(form, data = d, family = fam), error = function(e) e),
      warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
      message = function(c) invokeRestart("muffleMessage"))
    if (inherits(f, "error")) next
    # a plain fit of the same data with no latent bounds nothing; the same family
    # with no latent gives the scale of a sane logLik
    p <- tryCatch(as.numeric(logLik(suppressWarnings(frm(plain, data = d, family = fam)))), error = function(e) NA)
    par <- f$opt$par
    out[[length(out) + 1]] <- data.frame(seed = s, ll = as.numeric(logLik(f)), plain = p, code = f$opt$convergence, nw = length(w),
                                         ldisp = if (any(names(par) == "betad")) par[names(par) == "betad"][1] else NA)
  }
  o <- do.call(rbind, out)
  o$over <- o$ll - o$plain
  cat(sprintf("%-22s fits %d  code!=0 %d  ldisp>30 %d  median ll-plain %.2f  largest ll %.4g (seed %d, code %d, warnings %d, ldisp %.1f)  silent with ldisp>30: %d\n",
              lab, nrow(o), sum(o$code != 0), sum(o$ldisp > 30, na.rm = TRUE), median(o$over), max(o$ll), o$seed[which.max(o$ll)], o$code[which.max(o$ll)], o$nw[which.max(o$ll)], o$ldisp[which.max(o$ll)], sum(o$ldisp > 30 & o$nw == 0, na.rm = TRUE)))
  invisible(o)
}
run("beta_binomial, no od", function(d) rbinom(nrow(d), 10, plogis(.3 * d$x + d$e)), y | trials(tr) ~ x + ar(t, g, cov = TRUE), beta_binomial(), y | trials(tr) ~ x)
run("Gamma log, shape 1000", function(d) rgamma(nrow(d), 1000, 1000 / exp(.5 + .3 * d$x + d$e)), y ~ x + ar(t, g, cov = TRUE), Gamma(link = "log"), y ~ x)
run("Beta, phi 1000", function(d) { m <- plogis(.3 * d$x + d$e); pmin(pmax(rbeta(nrow(d), 1000 * m, 1000 * (1 - m)), 1e-4), 1 - 1e-4) }, y ~ x + ar(t, g, cov = TRUE), Beta(), y ~ x)
run("negbinomial, pois", function(d) rpois(nrow(d), exp(1 + .3 * d$x + d$e)), y ~ x + ar(t, g, cov = TRUE), negbinomial(), y ~ x)
run("weibull, shape 50", function(d) rweibull(nrow(d), 50, exp(.5 + .3 * d$x + d$e)), y ~ x + ar(t, g, cov = TRUE), weibull(), y ~ x)
run("zero_one_infl beta", function(d) { m <- plogis(.3 * d$x + d$e); v <- rbeta(nrow(d), 1000 * m, 1000 * (1 - m)); u <- runif(nrow(d)); ifelse(u < .1, rbinom(nrow(d), 1, .5), v) }, y ~ x + ar(t, g, cov = TRUE), zero_one_inflated_beta(), y ~ x)
