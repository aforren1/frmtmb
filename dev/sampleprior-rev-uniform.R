# Reviewer of lane sampleprior: uniform() on doubles and on the taped
# objectives (the sampler's, and a MAP fit's), inside and outside its
# support; its constant against -log(ub - lb).
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
core <- asNamespace("frmtmb")
u <- core$parse_prior_dist("uniform(-1, 2)")
cat("doubles: ", sapply(c(-1.5, -1, -0.5, 1.9, 2, 2.5), function(v)
  core$prior_base_logdens(v, u)), " target -log(3) =", -log(3), "\n")
# the sampler's tape, a uniform on the slope
d <- D$wide
pr <- set_prior("uniform(-1, 2)", class = "b", coef = "x")
od <- cap_obj(y ~ x, d, gaussian(), prior = pr)
o0 <- cap_obj(y ~ x, d, gaussian())        # the defaults alone
p <- od$env$par
at <- function(o, b) { q <- p; q[2] <- b; o$env$f(q) }
for (b in c(-1.5, -0.999, 0.5, 1.999, 2.5)) {
  cat(sprintf("sampler tape, b_x = %6.3f: nll %s, uniform part (nll diff) %s\n",
              b, format(at(od, b)), format(at(od, b) - at(o0, b))))
}
cat("gradient outside:", od$gr(replace(od$par, 2, 2.5)), "\n")
cat("gradient inside :", od$gr(replace(od$par, 2, 0.5)), " vs defaults only",
    o0$gr(replace(o0$par, 2, 0.5)), "\n")
# a MAP fit's own tape
ft <- suppressWarnings(frm(y ~ x, data = d, family = gaussian(),
                           prior = set_prior("uniform(-1, 4)", class = "b",
                                             coef = "x")))
q <- ft$obj$env$last.par.best
cat("MAP fit estimate b_x", q[2], "\n")
for (b in c(-2, 0.5, 4.5)) {
  qq <- q; qq[2] <- b
  cat(sprintf("MAP tape, b_x = %4.1f: fn %s\n", b, format(ft$obj$fn(qq))))
}
# a uniform as the only prior on a slope, sampled: does NUTS live with
# the NaN wall?
ds <- suppressMessages(suppressWarnings(frm_sample(
  y ~ x, data = d, family = gaussian(),
  prior = set_prior("uniform(2.2, 2.6)", class = "b", coef = "x"),
  chains = 2, iter = 1000, refresh = 0, seed = 5)))
v <- ds$draws[, "b_x"]
cat("sampled b_x under uniform(2.2, 2.6): range", range(v),
    " divergences", sum(rstan::get_divergent_iterations(ds$stanfit)), "\n")
