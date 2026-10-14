# Which bin-2P blocks dev/brmsport-bin2-gen.R ports through an adapter,
# and the adapter each needs. A 2P block not listed here is classified
# from the survey alone (dev/brmsport-bin2-survey.R).
#
# G3 (tests.log_lik.R) and G4 (tests.posterior_predict.R) go through the
# thin adapter of dev/brmsport-bin2-adapters.R: brms's internal
# log_lik_<fun>(i, prep) and posterior_predict_<fun>(i, prep) become
# brms_adapt_log_lik(fun, i, prep) and brms_adapt_posterior_predict(),
# which evaluate frmtmb's family of the same brms name at column i of
# the prep. A block is listed when its prep is univariate and its rows
# independent, which is all the adapter reads; the ordinal blocks read
# thresholds too, which each ordinal family object takes through its
# exported ord_lpdf_make() and ord_sim_make() (fix round 1, after the
# review's dev/remeasure-rev-ord.R). Left to the survey: the
# multivariate, ARMA, SAR and FCOR blocks (the row density is not the
# family's), the categorical blocks (category arrays), wiener
# (frmtmb.eam's) and custom families.
#
# G2 (tests.distributions.R): the two blocks that compare brms's
# unexported ordinal densities dcumulative(), dsratio(), dcratio() and
# dacat() with brms's own closed-form helpers go through
# brms_adapt_dord(), the same ordinal family objects. The link-map
# blocks (link_<ordinal>(), inv_link_<ordinal>() on arrays) have no
# frmtmb counterpart and stay with the survey.
#
# A family frmtmb does not have (discrete_weibull, negbinomial2, frechet,
# gen_extreme_value) inside a listed block is substituted too: the
# adapter's brmsfamily() call then refuses it in frmtmb's own words, and
# the row records that refusal. Porting the block whole keeps the
# binomial, poisson, gamma, weibull and other comparisons beside them,
# which leaving the block to the survey would lose.
#
# The reference side of every assertion stays brms's exported density or
# base R's. brms's text calls those unqualified, because brms's tests run
# inside its namespace, so the header binds the exported ones.

# brms's two inverse links, as brms 2.23.0 defines them
# (R/numeric-helpers.R), bound under the names its setup lines call; a
# copy of the body rather than plogis(), so the data both sides read are
# bit-identical to brms's
brms_adapt_pre_links <- c(
  "inv_logit <- function(x) 1 / (1 + exp(-x))",
  "inv_cloglog <- function(x) 1 - exp(-exp(x))")

brms_adapt_subst <- function(kind, funs) {
  a <- paste0("brms_adapt_", kind)
  c(
    # the direct call, brms:::log_lik_gaussian(1, prep = prep)
    lapply(funs, function(f) {
      c(paste0("brms:::", kind, "_", f, "("),
        paste0(a, "(\"", f, "\", "))
    }),
    # the sapply() form, sapply(1:nobs, brms:::log_lik_gaussian, ...)
    lapply(funs, function(f) {
      c(paste0("brms:::", kind, "_", f, ","),
        paste0(a, ", fun = \"", f, "\","))
    }),
    list(c("brms:::inv_logit(", "inv_logit("),
         c("brms:::inv_cloglog(", "inv_cloglog("),
         # brms exports rasym_laplace(); the internal spelling is the
         # same function
         c("brms:::rasym_laplace(", "rasym_laplace("),
         # brms skips on a package its OWN family code needs; frmtmb's
         # families need neither (dev/brmsport-bin2-adapters.R)
         c("skip_if_not_installed(\"extraDistr\")",
           "brms_adapt_brms_needs(\"extraDistr\")"),
         c("skip_if_not_installed(\"betareg\")",
           "brms_adapt_brms_needs(\"betareg\")")))
}

brms_adapt_funs <- c(
  "gaussian", "student", "lognormal", "shifted_lognormal", "exgaussian",
  "skew_normal", "asym_laplace", "zero_inflated_asym_laplace", "binomial",
  "beta_binomial", "discrete_weibull", "poisson", "negbinomial",
  "negbinomial2", "geometric", "com_poisson", "exponential", "gamma",
  "frechet", "inverse.gaussian", "weibull", "gen_extreme_value",
  "bernoulli", "beta", "von_mises", "hurdle_poisson", "hurdle_negbinomial",
  "hurdle_gamma", "zero_inflated_poisson", "zero_inflated_negbinomial",
  "zero_inflated_binomial", "zero_inflated_beta_binomial",
  "zero_inflated_beta", "zero_one_inflated_beta", "xbeta",
  "cumulative", "sratio", "cratio", "acat", "hurdle_cumulative",
  # fix round 2: the categorical blocks are tried with no new adapter
  # code; frmtmb's refusal is what a row that cannot run records
  "categorical", "multinomial", "dirichlet_multinomial", "dirichlet",
  "dirichlet2", "logistic_normal")

# brms's exported d, r and p functions that dev/brmsport-bin2-adapters.R
# shims (brms_g1_shims), bound by name in the distributions file's header
brms_g1_names <- c(
  "dstudent_t", "rstudent_t", "dvon_mises", "rvon_mises", "dskew_normal",
  "rskew_normal", "dexgaussian", "rexgaussian", "dinv_gaussian",
  "pinv_gaussian", "rinv_gaussian", "dbeta_binomial", "rbeta_binomial",
  "dasym_laplace", "rasym_laplace", "dzero_inflated_poisson",
  "dzero_inflated_negbinomial", "dzero_inflated_binomial",
  "dzero_inflated_beta_binomial", "dzero_inflated_beta",
  "dhurdle_poisson", "dhurdle_negbinomial", "dhurdle_gamma",
  "dhurdle_lognormal")

bin2_adapted <- list(
  "tests.log_lik.R" = list(
    topic = "log_lik", pkgs = "frmtmb",
    blocks = c(
      "log_lik for location shift models works as expected",
      "log_lik for various skewed normal models works as expected",
      "log_lik of aysm_laplace models runs without errors",
      "log_lik for count and survival models works correctly",
      "log_lik for bernoulli and beta models works correctly",
      "log_lik for circular models runs without errors",
      "log_lik for zero-inflated and hurdle models runs without erros",
      "censored and truncated log_lik run without errors",
      "log_lik for the xbeta model runs without errors",
      "log_lik for ordinal models runs without erros",
      "log_lik for categorical and related models runs without erros"),
    pre = c(
      "dstudent_t <- brms::dstudent_t",
      "dshifted_lnorm <- brms::dshifted_lnorm",
      "dexgaussian <- brms::dexgaussian",
      "dskew_normal <- brms::dskew_normal",
      "dbeta_binomial <- brms::dbeta_binomial",
      "dfrechet <- brms::dfrechet",
      "dinv_gaussian <- brms::dinv_gaussian",
      "dgen_extreme_value <- brms::dgen_extreme_value",
      "rasym_laplace <- brms::rasym_laplace",
      brms_adapt_pre_links,
      "SW <- suppressWarnings"),
    subst = brms_adapt_subst("log_lik", brms_adapt_funs)),
  "tests.posterior_predict.R" = list(
    topic = "posterior_predict", pkgs = "frmtmb",
    blocks = c(
      "posterior_predict for location shift models runs without errors",
      "posterior_predict for various skewed models runs without errors",
      "posterior_predict for aysm_laplace models runs without errors",
      paste("posterior_predict for count and survival models runs",
            "without errors"),
      "posterior_predict for bernoulli and beta models works correctly",
      "posterior_predict for xbeta models works correctly",
      "posterior_predict for circular models runs without errors",
      paste("posterior_predict for zero-inflated and hurdle models runs",
            "without erros"),
      "truncated posterior_predict run without errors",
      "posterior_predict for ordinal models runs without errors",
      paste("posterior_predict for categorical and related models runs",
            "without erros")),
    pre = brms_adapt_pre_links,
    subst = brms_adapt_subst("posterior_predict", brms_adapt_funs)),
  "tests.distributions.R" = list(
    topic = "distributions", pkgs = "frmtmb",
    blocks = c(
      "d<ordinal_family>() works correctly",
      paste("dsratio() and dcratio() give the same results for symmetric",
            "distribution functions"),
      # G1 (fix round 2): every d/p/q/r block but the multivariate normal
      # and t (no density outside a fit) and wiener (brms compares two of
      # its own backends)
      "student distribution works correctly",
      "von_mises distribution functions run without errors",
      "skew_normal distribution functions run without errors",
      "exgaussian distribution functions run without errors",
      "frechet distribution functions run without errors",
      "inv_gaussian distribution functions run without errors",
      "beta_binomial distribution functions run without errors",
      "gen_extreme_value distribution functions run without errors",
      "asym_laplace distribution functions run without errors",
      "zero-inflated distribution functions run without errors",
      "hurdle distribution functions run without errors",
      # G2 (fix 2b): the inv_link_<ordinal>() array maps
      "inv_link_<ordinal_family>() works correctly for arrays",
      paste("inv_link_sratio() and inv_link_cratio() applied to arrays",
            "give the same results for symmetric distribution functions",
            "(when respecting the sign appropriately).")),
    pre = c(sprintf("%s <- brms_g1_shims[[\"%s\"]]", brms_g1_names,
                    brms_g1_names),
            sprintf("inv_link_%s <- brms_g2_shims[[\"%s\"]]",
                    c("cumulative", "sratio", "cratio", "acat"),
                    c("cumulative", "sratio", "cratio", "acat")),
            "SW <- suppressWarnings"),
    subst = c(
      list(c("skip_if_not_installed(\"extraDistr\")",
             "brms_adapt_brms_needs(\"extraDistr\")")),
      lapply(c("cumulative", "sratio", "cratio", "acat"), function(f) {
        c(paste0("d", f, "(seq_len"),
          paste0("brms_adapt_dord(\"", f, "\", seq_len"))
      }),
      # testthat::test_path() would resolve against the staged run's
      # directory; identity() keeps the parentheses balanced
      list(c("testthat::test_path(file.path(\"helpers\", ",
             "identity(file.path(brms_adapt_suite_dir, \"helpers\", "))))
)
