# Reviewer of lane sampleprior, fix round 1: the guards' control on the
# SAMPLER's objective (the lane's control is on frm()'s): with 0.71.0's
# unguarded densities swapped in, the captured frm_sample() tape must be
# finite at b_x = 1.5 under beta(1, 1), and NaN with the guards; inside,
# value and gradient bitwise equal.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
core <- asNamespace("frmtmb")
guarded <- core$prior_base_logdens
unguarded <- function(x, dist) {
  switch(dist$kind,
    gamma = dist$shape * log(dist$rate) - lgamma(dist$shape) +
      (dist$shape - 1) * log(x) - dist$rate * x,
    inv_gamma = dist$shape * log(dist$scale) - lgamma(dist$shape) -
      (dist$shape + 1) * log(x) - dist$scale / x,
    beta = (dist$shape1 - 1) * log(x) + (dist$shape2 - 1) * log1p(-x) -
      (lgamma(dist$shape1) + lgamma(dist$shape2) -
         lgamma(dist$shape1 + dist$shape2)),
    guarded(x, dist))
}
grab <- function(f, ...) {
  utils::assignInNamespace("prior_base_logdens", f, "frmtmb")
  on.exit(utils::assignInNamespace("prior_base_logdens", guarded, "frmtmb"))
  cap_obj(...)
}
pr <- suppressWarnings(set_prior("beta(1, 1)", class = "b", coef = "x"))
o1 <- grab(guarded, y ~ x, D$unit, frmtmb::brmsfamily("beta"), prior = pr)
o0 <- grab(unguarded, y ~ x, D$unit, frmtmb::brmsfamily("beta"), prior = pr)
q <- o1$par; q[2] <- 1.5
cat("sampler tape, b_x = 1.5: guarded", o1$fn(q), " unguarded", o0$fn(q), "\n")
# the defaults of a zinb (zi beta(1, 1), shape inv_gamma) inside the support
o1 <- grab(guarded, y ~ x, D$zcnt, frmtmb::brmsfamily("zero_inflated_negbinomial"))
o0 <- grab(unguarded, y ~ x, D$zcnt, frmtmb::brmsfamily("zero_inflated_negbinomial"))
set.seed(2)
same <- vapply(1:20, function(k) {
  q <- o1$par + rnorm(length(o1$par), 0, 0.5)
  identical(o1$fn(q), o0$fn(q)) && identical(o1$gr(q), o0$gr(q))
}, NA)
cat("zinb defaults, 20 points: value and gradient bitwise equal at", sum(same), "of 20\n")
