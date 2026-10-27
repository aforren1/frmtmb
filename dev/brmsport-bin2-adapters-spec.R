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
# review's dev/remeasure-rev-ord.R). The categorical blocks (lane
# bin2port) map brms's per-category mu1, mu2, ... onto the categorical,
# multinomial, dirichlet, dirichlet2 and logistic_normal families
# (brms_adapt_cat_*()), and the wiener blocks onto frmtmb.eam's wiener
# family, with brms's dec as its dec() term. The SAR and FCOR blocks
# (lane loglik) go through brms_adapt_ac_log_lik() and
# brms_adapt_ac_posterior_predict(): brms's prep becomes a frmtmb model
# with one mean per row and the term over brms's matrix, and the values
# are frmtmb's autocor_row_loglik() and sim_draw(). ADAPTATION: the SAR
# blocks' Msar = diag(10) is replaced by the row-standardized path graph
# over the same 10 rows (brms_adapt_sar_path()). frmtmb refuses diag(10)
# by design: I - rho I makes rho and the scale one parameter, so sar()
# needs a negative and a positive eigenvalue. The value checks run
# brms's own function on the substituted matrix too. Left to the
# survey: the multivariate and ARMA blocks and custom families.
#
# G2 (tests.distributions.R): the two blocks that compare brms's
# unexported ordinal densities dcumulative(), dsratio(), dcratio() and
# dacat() with brms's own closed-form helpers go through
# brms_adapt_dord(), the same ordinal family objects; the
# inv_link_<ordinal>() array blocks through brms_g2_shims (fix 2b); and
# dcategorical() and inv_link_categorical() through brms_g2_cat_shims,
# frmtmb's categorical family (lane bin2port). The link maps the other
# way (link_<ordinal>(), link_categorical()) have no frmtmb counterpart
# and stay with the survey.
#
# A family frmtmb does not have (negbinomial2 and dirichlet_multinomial
# were two, until lane nb2dm) inside a listed block is substituted
# too: the adapter's brmsfamily()
# call then refuses it in frmtmb's own words, and the row records that
# refusal. Porting the block whole keeps the binomial, poisson, gamma,
# weibull and other comparisons beside them, which leaving the block to
# the survey would lose.
#
# Value checks (lane bin2port): where brms's assertion checks only a
# shape and the line before it computed the value through the adapter
# from a log density or a p function, the generator adds a check that
# compares that value with brms's own function on the same arguments
# (bin2_value_check_ok(), brms_value_check()). It is recorded as its
# own "value" row and never changes the brms row's verdict.
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
           "brms_adapt_brms_needs(\"betareg\")"),
         # lane bin2port: RWiener is brms's wiener backend, and brms
         # skips its wiener posterior_predict block for a fault of it
         c("skip_if_not_installed(\"RWiener\")",
           "brms_adapt_brms_needs(\"RWiener\")"),
         c("skip(\"skip as long as RWiener fails",
           "brms_adapt_brms_skip(\"skip as long as RWiener fails")))
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
  "dirichlet2", "logistic_normal",
  # lane bin2port: frmtmb.eam's wiener
  "wiener")

# brms's residual correlation helpers the SAR and FCOR blocks call
# (lane loglik), and the one adaptation they need: Msar = diag(10)
# becomes the path graph (see the header)
brms_adapt_ac_spec_funs <- c("gaussian_lagsar", "student_lagsar",
                             "gaussian_errorsar", "student_errorsar",
                             "gaussian_fcor", "student_fcor")
brms_adapt_sar_subst <- list(c("Msar = diag(10)",
                               "Msar = brms_adapt_sar_path(10)"))

# The rows that hold only on the substituted input. The ledger keeps
# them as passes, of class "adapted-input", with this reason, so that
# they are never counted as literal passes of brms's rows: on brms's own
# Msar = diag(10) frmtmb refuses the model, by design
# (dev/brmsport-bin2-ledger.R).
bin2_adapted_input_reason <- paste(
  "pass on an adapted input: brms's prep uses Msar = diag(10), and the",
  "port replaces it by the row-standardized path graph over the same 10",
  "rows (brms_adapt_sar_path()). frmtmb refuses brms's literal input by",
  "design: I - rho I makes rho and the scale one parameter, so sar()",
  "needs a weights matrix with a negative and a positive eigenvalue.",
  "brms's assertion (the dimensions) holds on the substituted matrix,",
  "and for log_lik the value check against brms's own function on it",
  "holds too.")
bin2_adapted_input <- stats::setNames(
  rep(bin2_adapted_input_reason, 8L),
  c("log_lik:150", "log_lik:152", "log_lik:157", "log_lik:159",
    "posterior_predict:111", "posterior_predict:113",
    "posterior_predict:118", "posterior_predict:120"))

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
  "dhurdle_lognormal",
  # lane bin2port: 0.71.0's extreme-value families
  "dfrechet", "pfrechet", "rfrechet", "dgen_extreme_value",
  "pgen_extreme_value", "rgen_extreme_value",
  # lane cdfs: the p functions of the families that gained a CDF
  "pstudent_t", "pvon_mises", "pskew_normal", "pexgaussian",
  "pbeta_binomial", "pasym_laplace", "pzero_inflated_poisson",
  "pzero_inflated_negbinomial", "pzero_inflated_binomial",
  "pzero_inflated_beta_binomial", "pzero_inflated_beta",
  "phurdle_poisson", "phurdle_negbinomial", "phurdle_gamma",
  "phurdle_lognormal")
brms_g2_cat_names <- c("dcategorical", "inv_link_categorical")

# What a pass asserts, read from brms's assertion text: a shape, length,
# names or type check, as against a comparison of values. The ledger
# classifies passes with it; the generator adds a value check after a
# shape assertion (dev/brmsport-bin2-ledger.R). The shape call must be
# the assertion's FIRST argument: expect_equal(rowSums(pred), rep(1,
# nrow(pred))) compares values, and a match anywhere counted it as a
# shape (review 2026-10-14, finding 5). SW() may wrap the assertion.
bin2_shape_re <- paste0(
  "^[[:space:]]*(SW[(][[:space:]]*)?(",
  "expect_(equal|true|identical|equivalent)[(][[:space:]]*",
  "(dim|length|names|dimnames|nrow|ncol|is[.]numeric|is[.]character)[(]|",
  "expect_is[(]|expect_output[(]|expect_length[(]|expect_named[(])")

# The value-check candidate at statement j of a block: the name the line
# before it assigns, when that line computes a log density or a CDF
# (brms:::log_lik_<family>(), or a d or p function the header binds to a
# shim) and statement j is one assertion of a shape of that name; NULL
# otherwise. A line that draws (r functions, posterior_predict_*) is not
# a candidate: two simulators agree in distribution, not draw by draw.
bin2_value_check_ok <- function(b, j) {
  if (n_expect(b$stmts[[j]]) != 1L || !grepl(bin2_shape_re, b$texts[j])) {
    return(NULL)
  }
  # the call inside the assertion itself, as in
  # expect_equal(length(brms:::log_lik_wiener(i, prep)), ns)
  inner <- bin2_inner_log_lik(b$stmts[[j]])
  if (!is.null(inner)) return(list(var = ".value", inline = inner))
  if (j < 2L || n_expect(b$stmts[[j - 1L]]) != 0L) return(NULL)
  s <- b$stmts[[j - 1L]]
  if (!is.call(s) || !as.character(s[[1]])[1] %in% c("<-", "=") ||
      !is.name(s[[2]])) {
    return(NULL)
  }
  v <- as.character(s[[2]])
  if (!v %in% all.vars(b$stmts[[j]])) return(NULL)
  rhs <- s[[3]]
  if (is.call(rhs) && identical(rhs[[1]], as.name("SW"))) rhs <- rhs[[2]]
  fn <- if (is.call(rhs)) deparse(rhs[[1]]) else ""
  dp <- grep("^[dp]", c(brms_g1_names, brms_g2_cat_names), value = TRUE)
  ref <- grepl("brms:::log_lik_", b$texts[j - 1L], fixed = TRUE) ||
    fn %in% dp
  if (ref) list(var = v, inline = NULL) else NULL
}

# The one brms:::log_lik_<family>() call inside an expression, or NULL
bin2_inner_log_lik <- function(e) {
  found <- list()
  walk <- function(x) {
    if (!is.call(x)) return(invisible())
    h <- deparse(x[[1]])
    if (startsWith(h, "brms:::log_lik_")) {
      found[[length(found) + 1L]] <<- x
      return(invisible())
    }
    # x[[k]] is read in place: an empty argument (x[i, ]) cannot be
    # bound to a name
    for (k in seq_along(x)[-1L]) if (is.call(x[[k]])) walk(x[[k]])
  }
  walk(e)
  if (length(found) == 1L) found[[1L]] else NULL
}

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
      "log_lik for categorical and related models runs without erros",
      # lane bin2port: frmtmb.eam's wiener family
      "log_lik for the wiener diffusion model runs without errors",
      # lane loglik: the residual correlation terms
      "log_lik for SAR models runs without errors",
      "log_lik for FCOR models runs without errors",
      "log_lik for FCOR student-t models works as expected"),
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
    subst = c(brms_adapt_subst("log_lik", c(brms_adapt_funs,
                                            brms_adapt_ac_spec_funs)),
              brms_adapt_sar_subst)),
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
            "without erros"),
      # lane bin2port: frmtmb.eam's wiener family
      paste("posterior_predict for the wiener diffusion model runs",
            "without errors"),
      # lane loglik: the residual correlation terms (brms's block label
      # for the SAR one says loglik)
      "loglik for SAR models runs without errors",
      "posterior_predict for FCOR models runs without errors"),
    pre = brms_adapt_pre_links,
    subst = c(brms_adapt_subst("posterior_predict",
                               c(brms_adapt_funs, brms_adapt_ac_spec_funs)),
              brms_adapt_sar_subst)),
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
            "(when respecting the sign appropriately)."),
      # lane bin2port: the categorical probability maps
      "dcategorical() works correctly",
      "inv_link_categorical() works correctly for arrays"),
    pre = c(sprintf("%s <- brms_g1_shims[[\"%s\"]]", brms_g1_names,
                    brms_g1_names),
            sprintf("inv_link_%s <- brms_g2_shims[[\"%s\"]]",
                    c("cumulative", "sratio", "cratio", "acat"),
                    c("cumulative", "sratio", "cratio", "acat")),
            sprintf("%s <- brms_g2_cat_shims[[\"%s\"]]", brms_g2_cat_names,
                    brms_g2_cat_names),
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
