# Short check of lane nlstop's fix round. Part 1: which families pass the
# scale gate (a betad entry in the outer vector), traced by a counter of
# optimizer runs inside prior_conflict_escape(). Part 2: the centered
# intercept fix against an intercept and a slope both in conflict, a dpar
# (sigma) intercept prior, 0 + Intercept, and an mvbf with a per-response
# intercept prior. Part 3: autoscale's seeded main fit.
# Usage: Rscript dev/reviews/nlstop-short-attack.R [part ...]
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
orig <- ns$prior_conflict_escape
cnt <- new.env(); cnt$k <- 0L; cnt$calls <- 0L; cnt$betad <- NA
oo <- ns$optimize_obj
utils::assignInNamespace("optimize_obj", function(...) {
  cnt$k <- cnt$k + 1L
  oo(...)
}, "frmtmb")
utils::assignInNamespace("prior_conflict_escape", function(obj, opt, ...) {
  cnt$calls <- cnt$calls + 1L
  cnt$betad <- any(names(obj$par) == "betad")
  k0 <- cnt$k
  o <- orig(obj, opt, ...)
  cnt$restarts <- cnt$k - k0
  o
}, "frmtmb")
reset <- function() {
  cnt$calls <- 0L; cnt$restarts <- NA; cnt$betad <- NA
}
show <- function(lab, r) {
  if (iserr(r)) { cat(sprintf("ROW %s ERR %s\n", lab, substr(r$msg, 1, 200))); return() }
  cat(sprintf("ROW %s llp %.6f code %s warn %s msg %s calls %d betad %s restarts %s\n",
              lab, llp(r), r$opt$convergence, wsum(r),
              any(attr(r, "mc") == "frmtmb_prior_multistart"), cnt$calls,
              cnt$betad, cnt$restarts))
  flush.console()
}
args <- commandArgs(TRUE)
want <- function(k) !length(args) || k %in% args

if (want("p1")) {
  set.seed(21); n <- 600
  d <- data.frame(x = rnorm(n), g = factor(sample(30, n, TRUE)))
  d$obs <- factor(seq_len(n))
  eta <- 0.5 + 2 * d$x
  pr <- set_prior("normal(0, 0.1)", class = "b")
  fams <- list(
    gaussian = list(y = 1 + 2 * d$x + rnorm(n), fam = gaussian()),
    student = list(y = 1 + 2 * d$x + rt(n, 4), fam = student()),
    negbinomial = list(y = rnbinom(n, mu = exp(eta), size = 3),
                       fam = negbinomial()),
    gamma = list(y = rgamma(n, 3, 3 / exp(eta)), fam = Gamma("log")),
    beta = list(y = plogis(eta / 4 + rnorm(n, 0, 0.5)), fam = Beta()),
    lognormal = list(y = exp(eta / 2 + rnorm(n, 0, 0.5)), fam = lognormal()),
    poisson = list(y = rpois(n, exp(eta)), fam = poisson()),
    bernoulli = list(y = rbinom(n, 1, plogis(eta)), fam = bernoulli()),
    zip = list(y = ifelse(runif(n) < 0.2, 0, rpois(n, exp(eta))),
               fam = zero_inflated_poisson()),
    cumulative = list(y = ordered(cut(eta + rlogis(n), c(-Inf, 0, 1, 2, Inf))),
                      fam = cumulative()))
  for (fn in names(fams)) {
    d$y <- fams[[fn]]$y
    reset(); show(fn, tt(frm(y ~ x, family = fams[[fn]]$fam, data = d, prior = pr)))
  }
  d$y <- 1 + 2 * d$x + rnorm(n, 0, exp(0.3 * d$x))
  reset(); show("gaussian_sigma_x", tt(frm(bf(y ~ x, sigma ~ x), data = d, prior = pr)))
  # a free scale only as a group-level sd: poisson with an observation-level
  # effect, slope 2 under normal(0, 0.1); default against a start with the
  # slope at 0 and the sd inflated
  d$y <- rpois(n, exp(eta + rnorm(n, 0, 0.3)))
  reset(); r0 <- tt(frm(y ~ x + (1 | obs), family = poisson(), data = d, prior = pr))
  show("poisson_olre", r0)
  r1 <- tt(frm(y ~ x + (1 | obs), family = poisson(), data = d, prior = pr,
               start = list(beta = c(0.5, 0), theta = log(2))))
  show("poisson_olre_start0", r1)
}

conflict_data <- function(seed) {
  set.seed(seed); n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$yl <- 1e6 * (1 + 0.05 * d$x) + rnorm(n, 0, 1e3)
  d
}
if (want("p2")) {
  d <- conflict_data(6201)
  # A: intercept and slope priors both in conflict
  prA <- c(set_prior("normal(0, 1e4)", class = "Intercept"),
           set_prior("normal(0, 1000)", class = "b"))
  reset(); rA <- tt(frm(yl ~ x, data = d, prior = prA)); show("A_default", rA)
  rA1 <- tt(frm(yl ~ x, data = d, prior = prA,
                start = list(beta = c(0, 0), betad = log(1e6))))
  show("A_start0", rA1)
  if (!iserr(rA)) cat("  A est:", format(fixef(rA)[, 1], digits = 7),
                      "sigma", exp(rA$estimates$betad[1]), "\n")
  # B: a prior on sigma's intercept (dpar), sigma's log near 6.9
  prB <- c(set_prior("normal(0, 0.5)", class = "Intercept", dpar = "sigma"),
           set_prior("normal(0, 1000)", class = "b"))
  reset(); rB <- tt(frm(bf(yl ~ x, sigma ~ 1), data = d, prior = prB))
  show("B_sigma_int", rB)
  # C: no centering
  prC <- set_prior("normal(0, 1e4)", class = "b", coef = "Intercept")
  reset(); rC <- tt(frm(yl ~ 0 + Intercept + x, data = d, prior = prC))
  show("C_0+Intercept", rC)
  rC1 <- tt(frm(yl ~ 0 + Intercept + x, data = d, prior = prC,
                start = list(beta = c(0, 5e4), betad = log(1e6))))
  show("C_start0", rC1)
  # D: mvbf, per-response intercept prior on yl only; no residual
  # correlation, so its penalized value is the sum of the two univariate ones
  set.seed(7); d$y2 <- 1e6 * (1 + 0.05 * d$x) + rnorm(nrow(d), 0, 1e3)
  prD <- set_prior("normal(0, 1e4)", class = "Intercept", resp = "yl")
  reset(); rD <- tt(frm(mvbf(bf(yl ~ x), bf(y2 ~ x), rescor = FALSE), data = d,
                        prior = prD))
  show("D_mvbf", rD)
  u1 <- tt(frm(yl ~ x, data = d,
               prior = set_prior("normal(0, 1e4)", class = "Intercept")))
  u2 <- tt(frm(y2 ~ x, data = d))
  cat(sprintf("  D univariate sum %.6f (yl %.6f, y2 %.6f)\n", llp(u1) + llp(u2),
              llp(u1), llp(u2)))
  if (!iserr(rD)) print(fixef(rD)[, 1])
}
if (want("p3")) {
  # autoscale: x on a 1e-4 scale qualifies; slope 5e8 under normal(0, 1e7)
  d <- conflict_data(6201)
  d$xs <- d$x * 1e-4
  pr <- c(set_prior("normal(0, 1e7)", class = "b"),
          set_prior("normal(1e6, 1e6)", class = "Intercept"))
  for (as in list(NULL, TRUE, FALSE)) {
    reset()
    r <- tt(frm(yl ~ xs, data = d, prior = pr,
                control = frmtmb_control(autoscale = as)))
    show(paste0("autoscale_", if (is.null(as)) "NULL" else as), r)
  }
}
cat("DONE\n")
