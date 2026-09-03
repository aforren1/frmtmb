# Translation of brms's "Parameterization of Response Distributions in
# brms" (brms 2.23.0, doc/brms_families.Rmd) onto the frmtmb surface.
#
# This vignette has no runnable model chunks: it is a reference that
# states, in mathematics, what each family's density is and which
# parameter the linear predictor drives. So the translation is a
# COVERAGE script. For every family the vignette names, one row records
# three things:
#
#   1. does frmtmb accept the constructor under brms's own name,
#   2. does its default link match brms's documented default, and
#   3. does the density frmtmb actually evaluates carry the
#      parameterization the vignette states.
#
# Point 3 is measured, not read. Each check integrates or sums the
# family's OWN lpdf and compares the result against the vignette's
# stated mean and standard deviation, so a family that documented one
# parameterization and computed another would fail here.
#
# Every label carries the "ML: " prefix. There is no sampling section:
# with no models in the vignette there is nothing to sample. One block
# at the end records what frm_sample() does with a representative
# family from each shape.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_families")

options(mc.cores = 1)
set.seed(1234)

## ---- the checkers ---------------------------------------------------

# frmtmb's by-name family lookup. It is internal because a user reaches
# it through frm(family =), but the audit needs the resolved object
# without paying for a fit.
FAM <- function(nm) frmtmb:::as_frmtmb_family(nm)

# The default link of the frmtmb constructor, against brms's own.
# `brms::brmsfamily(nm)$link` is the documented default in executable
# form, so the comparison is measured rather than transcribed.
link_ok <- function(nm, got) {
  want <- brms::brmsfamily(nm)$link
  if (!identical(got, want)) {
    stop(sprintf("link default differs: frmtmb '%s', brms '%s'", got, want))
  }
  invisible(TRUE)
}

# The mu link of a resolved family. Ordinal and categorical families
# apply their link inside the density rather than to mu, so their
# constructor's own argument default is the honest thing to compare.
mu_link <- function(fam) fam$links[[1]]$name
arg_link <- function(ctor) eval(formals(ctor)$link)

dens <- function(fam, y, dpars, aterms = list()) {
  exp(fam$lpdf(y, dpars, aterms))
}

# A continuous family: the density must integrate to one over its
# support, and its first two moments must be the ones the vignette
# names.
cont_ok <- function(fam, dpars, lo, hi, mean_want, sd_want = NULL,
                    aterms = list(), tol = 1e-4) {
  f <- function(y) dens(fam, y, dpars, aterms)
  I <- function(g) stats::integrate(g, lo, hi, subdivisions = 4000L,
                                    rel.tol = 1e-9)$value
  tot <- I(f)
  m <- I(function(y) y * f(y))
  stopifnot(abs(tot - 1) < 1e-5, abs(m - mean_want) < tol * max(1, abs(mean_want)))
  if (!is.null(sd_want)) {
    v <- I(function(y) y^2 * f(y)) - m^2
    stopifnot(abs(sqrt(v) - sd_want) < tol * max(1, sd_want))
  }
  sprintf("integral %.6f, mean %.5f (want %.5f)", tot, m, mean_want)
}

# A discrete family: the same two checks, summed over the support.
disc_ok <- function(fam, dpars, ys, mean_want, sd_want = NULL,
                    aterms = list(), tol = 1e-4) {
  p <- dens(fam, ys, dpars, aterms)
  tot <- sum(p)
  m <- sum(ys * p)
  stopifnot(abs(tot - 1) < 1e-6, abs(m - mean_want) < tol * max(1, abs(mean_want)))
  if (!is.null(sd_want)) {
    v <- sum(ys^2 * p) - m^2
    stopifnot(abs(sqrt(v) - sd_want) < tol * max(1, sd_want))
  }
  sprintf("total %.8f, mean %.5f (want %.5f)", tot, m, mean_want)
}

# The families the vignette names but frmtmb does not have. Each one
# still gets a row, so the count is complete.
absent <- function(nm) FAM(nm)

## ============ PATH 1: ML / Laplace (frm) ============

## ---- Location shift models ------------------------------------------

bv("post", "ML: gaussian", {
  link_ok("gaussian", mu_link(FAM("gaussian")))
  cont_ok(FAM("gaussian"), list(mu = 2, sigma = 3), -60, 70, 2, 3)
}, NA_character_, "")

bv("post", "ML: student", {
  link_ok("student", mu_link(FAM("student")))
  # nu = 6 so the variance exists: sd = sigma * sqrt(nu / (nu - 2)).
  cont_ok(FAM("student"), list(mu = 2, sigma = 3, nu = 6), -200, 200,
          2, 3 * sqrt(6 / 4))
}, NA_character_, "")

bv("post", "ML: skew_normal", {
  link_ok("skew_normal", mu_link(FAM("skew_normal")))
  # The vignette's whole point about skew_normal is that omega and xi
  # are reparameterized so mu is the mean and sigma the sd. Both hold.
  cont_ok(FAM("skew_normal"), list(mu = 2, sigma = 3, alpha = 5),
          -60, 70, 2, 3)
}, NA_character_, "")

## ---- Binary and count data models -----------------------------------

bv("post", "ML: binomial", {
  link_ok("binomial", mu_link(FAM("binomial")))
  disc_ok(FAM("binomial"), list(mu = 0.3), 0:10, 10 * 0.3,
          sqrt(10 * 0.3 * 0.7), aterms = list(trials = 10))
}, NA_character_, "")

bv("post", "ML: bernoulli", {
  link_ok("bernoulli", mu_link(bernoulli()))
  disc_ok(bernoulli(), list(mu = 0.3), 0:1, 0.3, sqrt(0.3 * 0.7))
}, NA_character_, "")

bv("post", "ML: poisson", {
  link_ok("poisson", mu_link(FAM("poisson")))
  disc_ok(FAM("poisson"), list(mu = 4), 0:80, 4, 2)
}, NA_character_, "")

bv("post", "ML: negbinomial", {
  link_ok("negbinomial", mu_link(negbinomial()))
  # The vignette's phi is frmtmb's `shape`: var = mu + mu^2 / phi.
  disc_ok(negbinomial(), list(mu = 4, shape = 2), 0:400, 4,
          sqrt(4 + 16 / 2))
}, "SPELLING", "the vignette's precision parameter phi is named shape in frmtmb, as in brms; the density is the same")

bv("post", "ML: geometric", {
  link_ok("geometric", mu_link(geometric()))
  # geometric is negbinomial with phi fixed to 1, exactly as stated.
  disc_ok(geometric(), list(mu = 4), 0:2000, 4, sqrt(4 + 16))
}, NA_character_, "")

# The vignette's discrete_weibull section is inside an HTML comment, so
# it is documented but not rendered. It is counted anyway: it is a name
# the vignette's source carries and brms implements.
bv("post", "ML: discrete_weibull", absent("discrete_weibull"),
   "MISSING", "no frmtmb path and no near substitute; the refusal lists every supported family, which is the right pointer")

## ---- Time-to-event models -------------------------------------------

bv("post", "ML: lognormal", {
  link_ok("lognormal", mu_link(lognormal()))
  # mu is the location on the LOG scale, so the response mean is
  # exp(mu + sigma^2 / 2). That is what the vignette's density says.
  cont_ok(lognormal(), list(mu = 1, sigma = 0.5), 1e-10, 200,
          exp(1 + 0.125), tol = 1e-3)
}, NA_character_, "")

# The one genuine link-default divergence among the families frmtmb
# has. brms documents Gamma with a log link; frmtmb's by-name route
# agrees, but stats::Gamma() carries link = "inverse" and frmtmb takes
# it at face value, exactly as brms's own validate_family() does. So
# `family = "Gamma"` and `family = Gamma()` are two different models in
# both packages, which is a shared trap rather than a frmtmb one.
bv("post", "ML: Gamma", {
  link_ok("Gamma", mu_link(FAM("Gamma")))
  stopifnot(identical(mu_link(frmtmb:::as_frmtmb_family(stats::Gamma())),
                      "inverse"))
  cont_ok(FAM("Gamma"), list(mu = 3, shape = 4), 1e-10, 120, 3, 3 / 2)
}, "BEHAVIOR",
"family = 'Gamma' gives the log link brms documents, but family = Gamma() gives inverse, because the stats constructor's own default is carried through; brms behaves the same way, so a reader is caught by it in both")

bv("post", "ML: weibull", {
  link_ok("weibull", mu_link(weibull()))
  # The vignette's s = mu / Gamma(1 + 1/alpha) reparameterization, so
  # mu is the mean, is present verbatim in the lpdf.
  cont_ok(weibull(), list(mu = 3, shape = 2), 1e-10, 60, 3,
          3 * sqrt(gamma(2) / gamma(1.5)^2 - 1))
}, NA_character_, "")

bv("post", "ML: exponential", {
  link_ok("exponential", mu_link(exponential()))
  cont_ok(exponential(), list(mu = 3), 1e-10, 400, 3, 3)
}, NA_character_, "")

# brms documents inverse.gaussian with the 1/mu^2 link, and stats
# supplies exactly that. frmtmb has no such link, so the ordinary
# spelling is refused; the by-name route silently gives log instead.
bv("post", "ML: inverse.gaussian", {
  link_ok("inverse.gaussian", mu_link(FAM("inverse.gaussian")))
}, "BEHAVIOR",
"frmtmb's inverse.gaussian defaults to log where brms defaults to 1/mu^2, and stats::inverse.gaussian() is refused outright with \"Unknown link: '1/mu^2'\", which names the available links but not that a different default is in force")

bv("post", "ML: inverse.gaussian mean check under the log link", {
  cont_ok(inverse.gaussian_fam <- FAM("inverse.gaussian"),
          list(mu = 3, shape = 5), 1e-8, 400, 3, tol = 1e-3)
}, NA_character_, "")

# cox is present now. (The v0.34 audit recorded it as missing.)
bv("post", "ML: cox", {
  link_ok("cox", mu_link(cox()))
  # The M-spline baseline hazard means the density needs the spline
  # basis columns that only a fit builds, so the check here is the
  # link and the presence of cox_baseline(), not a numeric integral.
  stopifnot(is.function(cox_baseline))
  "constructor and link present; baseline hazard is fit-time"
}, "BEHAVIOR", "the proportional-hazards form matches, but the baseline hazard is built at fit time from cox_baseline(), so it cannot be checked from the family object alone")

## ---- Extreme value models -------------------------------------------

bv("post", "ML: frechet", absent("frechet"), "MISSING",
   "no frmtmb path; weibull() is the nearest relative the vignette itself names but it is a different distribution")

bv("post", "ML: gen_extreme_value", absent("gen_extreme_value"), "MISSING",
   "no frmtmb path and no substitute; the refusal lists the supported families and none of them generalize weibull and frechet")

## ---- Response time models -------------------------------------------

bv("post", "ML: exgaussian", {
  link_ok("exgaussian", mu_link(exgaussian()))
  # The vignette's mu = xi + beta reparameterization is in the lpdf:
  # dexgauss(mu - beta, sigma, 1/beta). So mu is the mean.
  cont_ok(exgaussian(), list(mu = 5, sigma = 1, beta = 2), -20, 120, 5,
          sqrt(1 + 4))
}, NA_character_, "")

bv("post", "ML: shifted_lognormal", {
  link_ok("shifted_lognormal", mu_link(shifted_lognormal()))
  # lognormal shifted right by ndt, exactly as stated.
  cont_ok(shifted_lognormal(), list(mu = 1, sigma = 0.5, ndt = 0.3),
          0.3 + 1e-10, 200, 0.3 + exp(1 + 0.125), tol = 1e-3)
}, NA_character_, "")

bv("post", "ML: wiener", absent("wiener"), "MISSING",
   "no frmtmb path; the four-parameter diffusion model with its bs/ndt/bias dpars has no analog, and none is named")

## ---- Quantile regression --------------------------------------------

bv("post", "ML: asym_laplace", {
  link_ok("asym_laplace", mu_link(asym_laplace()))
  fam <- asym_laplace()
  dp <- list(mu = 2, sigma = 1.5, quantile = 0.8)
  f <- function(y) dens(fam, y, dp)
  # The defining property the vignette states: P(Y < mu) = p.
  below <- stats::integrate(f, -200, 2, subdivisions = 4000L)$value
  stopifnot(abs(below - 0.8) < 1e-6)
  sprintf("P(Y < mu) = %.8f (want 0.8)", below)
}, "SPELLING", "the vignette's quantile parameter p is the dpar named `quantile`, as in brms")

## ---- Probability models ---------------------------------------------

bv("post", "ML: Beta", {
  link_ok("beta", mu_link(Beta()))
  # var = mu (1 - mu) / (1 + phi), which is the vignette's density.
  cont_ok(Beta(), list(mu = 0.3, phi = 8), 1e-10, 1 - 1e-10, 0.3,
          sqrt(0.3 * 0.7 / 9))
}, "SPELLING", "the constructor is Beta() with a capital B, as in brms, but the family object names itself 'beta' in lower case, so summary() prints a name the constructor does not have")

bv("post", "ML: dirichlet", absent("dirichlet"), "MISSING",
   "no frmtmb path for a simplex-valued response; multinomial() takes counts, not proportions, so it is not a substitute and nothing names one")

bv("post", "ML: logistic_normal", absent("logistic_normal"), "MISSING",
   "no frmtmb path; the alternative to dirichlet is missing for the same reason")

## ---- Circular models ------------------------------------------------

# von_mises is present now. (The v0.34 audit recorded it as missing.)
bv("post", "ML: von_mises", {
  link_ok("von_mises", mu_link(von_mises()))
  fam <- von_mises()
  dp <- list(mu = 0.5, kappa = 3)
  f <- function(y) dens(fam, y, dp)
  tot <- stats::integrate(f, -pi, pi, subdivisions = 4000L)$value
  # The circular mean, which is what mu means on a circle.
  cm <- atan2(stats::integrate(function(y) sin(y) * f(y), -pi, pi)$value,
              stats::integrate(function(y) cos(y) * f(y), -pi, pi)$value)
  stopifnot(abs(tot - 1) < 1e-6, abs(cm - 0.5) < 1e-6)
  sprintf("integral %.8f, circular mean %.6f (want 0.5)", tot, cm)
}, NA_character_, "")

## ---- Ordinal and categorical models ---------------------------------

# The four ordinal families apply their link inside the density, so
# their family object reports links$mu = "identity". The constructor's
# own `link` argument is the documented default and it matches brms.
for (nm in c("cumulative", "sratio", "cratio", "acat")) {
  local({
    n <- nm
    bv("post", paste0("ML: ", n), {
      ctor <- get(n, envir = as.environment("package:frmtmb"))
      link_ok(n, arg_link(ctor))
      fam <- ctor()
      # The vignette's density for cumulative is
      # g(tau_{y+1} - eta) - g(tau_y - eta). Checking that the K
      # category probabilities sum to one at a fixed threshold vector
      # is the strongest claim available without a fit.
      extra <- list(tau_raw = c(-1, log(1), log(1)))
      p <- exp(fam$lpdf(1:4, list(mu = 0.4), list(), extra))
      stopifnot(abs(sum(p) - 1) < 1e-8)
      sprintf("K = 4 probabilities sum to %.10f", sum(p))
    }, NA_character_, "")
  })
}

# categorical() is present now. (The v0.34 audit recorded it as
# missing, with multinomial() plus a matrix response as the substitute.)
bv("post", "ML: categorical", {
  link_ok("categorical", arg_link(categorical))
  fam <- categorical(K = 3)
  p <- exp(fam$lpdf(1:3, list(mu2 = 0.4, mu3 = -0.2), list()))
  stopifnot(abs(sum(p) - 1) < 1e-8)
  sprintf("K = 3 probabilities sum to %.10f, reference category is 1", sum(p))
}, NA_character_, "")

# multinomial is the one family that will not construct on brms's
# spelling: brms reads K off the response matrix and frmtmb demands it
# up front.
bv("post", "ML: multinomial [brms spelling]", multinomial(), "SPELLING",
   "multinomial() needs K; the refusal says so and shows multinomial(K = 3), so it points right")

bv("post", "ML: multinomial(K = 3)", {
  fam <- multinomial(K = 3)
  stopifnot(identical(fam$primary_dpars, c("mu2", "mu3")))
  paste("dpars:", paste(fam$dpars, collapse = ", "))
}, NA_character_, "")

## ---- Zero-inflated and hurdle models --------------------------------

bv("post", "ML: zero_inflated_poisson", {
  link_ok("zero_inflated_poisson", mu_link(zero_inflated_poisson()))
  # The vignette: f_z(0) = z + (1-z) f(0), f_z(y) = (1-z) f(y).
  # So the mean is (1 - z) mu.
  disc_ok(zero_inflated_poisson(), list(mu = 4, zi = 0.3), 0:80,
          0.7 * 4)
}, NA_character_, "")

bv("post", "ML: zero_inflated_binomial", {
  link_ok("zero_inflated_binomial", mu_link(zero_inflated_binomial()))
  disc_ok(zero_inflated_binomial(), list(mu = 0.3, zi = 0.2), 0:10,
          0.8 * 10 * 0.3, aterms = list(trials = 10))
}, NA_character_, "")

bv("post", "ML: zero_inflated_negbinomial", {
  link_ok("zero_inflated_negbinomial",
          mu_link(zero_inflated_negbinomial()))
  disc_ok(zero_inflated_negbinomial(), list(mu = 4, shape = 2, zi = 0.3),
          0:600, 0.7 * 4)
}, NA_character_, "")

bv("post", "ML: zero_inflated_beta", {
  link_ok("zero_inflated_beta", mu_link(zero_inflated_beta()))
  fam <- zero_inflated_beta()
  dp <- list(mu = 0.3, phi = 8, zi = 0.2)
  # The zero is an atom, so the check splits: the atom's mass plus the
  # continuous part must be one.
  atom <- exp(fam$lpdf(0, dp, list()))
  cont <- stats::integrate(function(y) dens(fam, y, dp), 1e-10,
                           1 - 1e-10, subdivisions = 4000L)$value
  stopifnot(abs(atom - 0.2) < 1e-8, abs(atom + cont - 1) < 1e-6)
  sprintf("atom %.6f + continuous %.6f = %.8f", atom, cont, atom + cont)
}, NA_character_, "")

bv("post", "ML: hurdle_poisson", {
  link_ok("hurdle_poisson", mu_link(hurdle_poisson()))
  # The vignette: f_z(0) = z, f_z(y) = (1-z) f(y) / (1 - f(0)).
  fam <- hurdle_poisson()
  dp <- list(mu = 4, hu = 0.3)
  p <- exp(fam$lpdf(0:80, dp, list()))
  stopifnot(abs(p[1] - 0.3) < 1e-10, abs(sum(p) - 1) < 1e-8)
  sprintf("P(0) = %.8f (want 0.3), total %.10f", p[1], sum(p))
}, NA_character_, "")

bv("post", "ML: hurdle_negbinomial", absent("hurdle_negbinomial"),
   "MISSING", "no frmtmb path, and it is the cheapest of the missing families: hurdle_poisson() and negbinomial() are both present and it is their composition")

bv("post", "ML: hurdle_gamma", {
  link_ok("hurdle_gamma", mu_link(hurdle_gamma()))
  fam <- hurdle_gamma()
  dp <- list(mu = 3, shape = 4, hu = 0.3)
  atom <- exp(fam$lpdf(0, dp, list()))
  cont <- stats::integrate(function(y) dens(fam, y, dp), 1e-10, 200,
                           subdivisions = 4000L)$value
  stopifnot(abs(atom - 0.3) < 1e-10, abs(atom + cont - 1) < 1e-6)
  sprintf("atom %.6f + continuous %.6f = %.8f", atom, cont, atom + cont)
}, NA_character_, "")

bv("post", "ML: hurdle_lognormal", {
  link_ok("hurdle_lognormal", mu_link(hurdle_lognormal()))
  fam <- hurdle_lognormal()
  dp <- list(mu = 1, sigma = 0.5, hu = 0.3)
  atom <- exp(fam$lpdf(0, dp, list()))
  cont <- stats::integrate(function(y) dens(fam, y, dp), 1e-10, 300,
                           subdivisions = 4000L)$value
  stopifnot(abs(atom - 0.3) < 1e-10, abs(atom + cont - 1) < 1e-6)
  sprintf("atom %.6f + continuous %.6f = %.8f", atom, cont, atom + cont)
}, NA_character_, "")

bv("post", "ML: zero_one_inflated_beta", absent("zero_one_inflated_beta"),
   "MISSING", "no frmtmb path, and the second cheap one: zero_inflated_beta() is present and this adds one more atom at y = 1")

## ---- The coverage tally ---------------------------------------------

bv("post", "ML: family coverage tally", {
  named <- c("gaussian", "student", "skew_normal", "binomial", "bernoulli",
             "poisson", "negbinomial", "geometric", "discrete_weibull",
             "lognormal", "Gamma", "weibull", "exponential",
             "inverse.gaussian", "cox", "frechet", "gen_extreme_value",
             "exgaussian", "shifted_lognormal", "wiener", "asym_laplace",
             "Beta", "dirichlet", "logistic_normal", "von_mises",
             "cumulative", "sratio", "cratio", "acat", "categorical",
             "multinomial", "zero_inflated_poisson",
             "zero_inflated_binomial", "zero_inflated_negbinomial",
             "zero_inflated_beta", "hurdle_poisson", "hurdle_negbinomial",
             "hurdle_gamma", "hurdle_lognormal", "zero_one_inflated_beta")
  # multinomial resolves only with K, so it is counted separately: the
  # family is present and the call is not brms's.
  respell <- "multinomial"
  ok <- vapply(named, function(n) {
    !inherits(try(FAM(if (identical(n, "Beta")) "beta" else n),
                  silent = TRUE), "try-error")
  }, TRUE)
  gone <- setdiff(names(ok)[!ok], respell)
  cat("named by the vignette:    ", length(named), "\n")
  cat("accepted on brms's call:  ", sum(ok), "\n")
  cat("accepted, different call: ", length(respell), " ", respell, "\n")
  cat("missing:                  ", length(gone), "\n")
  print(gone)
  cat("link default differs from brms: inverse.gaussian\n")
  cat("parameterization differs from brms: none\n")
  sum(ok)
}, NA_character_, "")

## ---- What frm_sample() makes of them --------------------------------
#
# There is no sampling section because there are no models to sample.
# This one block asks the narrower question the maintainer wants on the
# record: does the sampler accept a family that frm() accepts. Three
# representatives cover the three shapes the vignette organizes by,
# rather than all 32, because the answer is a property of the tape and
# not of the density.

bv("model", "ML: frm_sample accepts a continuous, a count and an ordinal family", {
  d <- data.frame(x = rnorm(60))
  d$y <- exp(1 + 0.5 * d$x + rnorm(60, sd = 0.3))
  d$k <- rpois(60, exp(1 + 0.5 * d$x))
  d$o <- factor(as.integer(cut(d$x, 3)), ordered = TRUE)
  out <- character(0)
  for (fam in list(lognormal(), poisson(), cumulative())) {
    resp <- switch(fam$family, lognormal = "y", poisson = "k", "o")
    f <- frm(stats::as.formula(paste(resp, "~ x")), data = d, family = fam)
    s <- frm_sample(f, chains = 1, iter = 300, warmup = 150,
                    seed = 1, cores = 1, refresh = 0)
    out <- c(out, sprintf("%s: %d draws", fam$family, ndraws(s)))
  }
  print(out)
}, NA_character_, "")

bv_done()
