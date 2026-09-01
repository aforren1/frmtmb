# Feature compatibility registry.
#
# The package's worst historical defects were pairwise feature
# combinations that were accepted and then quietly ignored (rescor
# dropping cens(), profile misreading bounds). A guard that does not
# exist looks exactly like a guard that passed, so "no error" was never
# evidence of support. This registry makes the third state explicit:
# every pair is declared to work, declared to be refused, or declared
# untested. Nothing is silently assumed to work.
#
# The registry is declarative data, not the guards themselves. Wiring
# the guards to read from it is a later migration. Until then the tests
# in test-compat.R check the declarations against the real behavior for
# a sample of pairs, and check that every named feature still exists.

# ---------------------------------------------------------------- features

# name: how the feature is written in a formula or a call argument.
# key:  the identifier the package itself uses, so tests can look the
#       feature up in family_registry, covstruct_registry, and the
#       parser vocabularies. Display names carry "()" for callable
#       features; the key never does.
frmtmb_compat_features_tbl <- function() {
  f <- function(name, kind) {
    key <- sub("\\(\\)$", "", name)
    if (key %in% names(frmtmb_compat_special_keys)) {
      key <- unname(frmtmb_compat_special_keys[[key]])
    }
    data.frame(name = name, key = key, kind = kind,
               stringsAsFactors = FALSE)
  }
  fams <- c("gaussian", "student", "lognormal", "shifted_lognormal",
            "skew_normal", "exgaussian", "asym_laplace", "Gamma",
            "weibull", "exponential", "inverse.gaussian", "beta",
            "tweedie", "poisson", "negbinomial", "nbinom1", "geometric",
            "compois", "binomial", "bernoulli", "beta_binomial",
            "multinomial", "zero_inflated_poisson",
            "zero_inflated_negbinomial", "zero_inflated_binomial",
            "zero_inflated_beta", "hurdle_poisson", "hurdle_gamma",
            "hurdle_lognormal", "cumulative", "sratio", "cratio",
            "acat")
  covs <- c("us", "diag", "homdiag", "cs", "ar1", "hetar1", "ou",
            "toep", "homtoep", "homcs", "exp", "gau", "mat", "rr",
            "equalto", "gr_cov", "gr_prec", "smooth", "gp", "hsgp")
  do.call(rbind, c(
    lapply(fams, f, kind = "family"),
    lapply(covs, f, kind = "covstruct"),
    lapply(c("weights()", "trials()", "cens()", "trunc()", "se()",
             "mi()", "vint()", "vreal()"), f, kind = "aterm"),
    lapply(c("s()", "t2()", "mo()", "mi_pred()", "gp_pred()",
             "cs_pred()"), f, kind = "special"),
    lapply(c("REML", "quadrature", "profile", "autoscale", "sparse_x",
             "priors", "bounds", "verbose"), f, kind = "mode"),
    lapply(c("mvbf", "rescor", "|ID|", "nl", "mixture",
             "mixture_mvn"), f, kind = "structure"),
    lapply(c("fitted", "predict", "simulate", "residuals",
             "residuals_osa", "emmeans", "frm_sample",
             "confint_profile", "hypothesis_profile"), f,
           kind = "method"),
    # formula-grammar spellings, which have their own restrictions and
    # belong in the table even though they name no package object
    lapply(c("bar_crossing", "call_group", "double_bar"), f,
           kind = "grammar")
  ))
}

# The special vocabulary collides with the covariance vocabulary: gp,
# cs and mi name both a bar-term structure and a predictor special.
# Display names disambiguate; these keys map back to the parser.
frmtmb_compat_special_keys <- c(mi_pred = "mi", gp_pred = "gp",
                                cs_pred = "cs")

# ------------------------------------------------------------------ groups

# Named feature sets, so one rule can cover "the families with a CDF"
# instead of thirty near-identical rows.
frmtmb_compat_groups_lst <- list(
  # the six families carrying an AD log-CDF, which is what cens() and
  # trunc() need to form their likelihood contributions
  cdf = c("gaussian", "poisson", "lognormal", "exponential", "weibull",
          "inverse.gaussian"),
  # cens() additionally refuses discrete responses, so poisson drops out
  cdf_continuous = c("gaussian", "lognormal", "exponential", "weibull",
                     "inverse.gaussian"),
  gaussian_like = c("gaussian", "student"),
  ordinal = c("cumulative", "sratio", "cratio", "acat"),
  # cs() category-specific effects are undefined under the cumulative
  # parameterization; the sequential and adjacent-category ones take them
  ordinal_cs = c("sratio", "cratio", "acat"),
  discrete = c("poisson", "negbinomial", "nbinom1", "geometric",
               "compois", "binomial", "bernoulli", "beta_binomial",
               "zero_inflated_poisson", "zero_inflated_negbinomial",
               "zero_inflated_binomial", "hurdle_poisson"),
  no_simulator = c("tweedie", "compois", "hurdle_poisson", "cumulative",
                   "sratio", "cratio", "acat", "multinomial"),
  matrix_response = c("multinomial"),
  trials_families = c("binomial", "beta_binomial",
                      "zero_inflated_binomial", "multinomial"),
  # quadrature marginalizes one scalar random effect at a time
  quadrature_blocks = c("us", "diag", "homdiag"),
  # every covariance structure quadrature cannot marginalize
  wide_blocks = c("cs", "ar1", "hetar1", "ou", "toep", "homtoep",
                  "homcs", "exp", "gau", "mat", "rr", "equalto",
                  "gr_cov", "gr_prec", "smooth", "gp", "hsgp"),
  # structures over a metric coordinate rather than positional levels
  spatial = c("exp", "gau", "mat"),
  # positional structures: the level order sets the lag, not the value
  positional = c("ar1", "hetar1", "toep", "homtoep"),
  latent_gp = c("gp", "hsgp"),
  known_cov = c("gr_cov", "gr_prec", "equalto"),
  post_fit = c("fitted", "predict", "simulate", "residuals",
               "residuals_osa", "emmeans")
)

# ------------------------------------------------------------------- rules

# Patterns, from least to most specific:
#   "*"                any feature
#   "kind:<kind>"      any feature of that kind
#   "group:<group>"    any feature in that named group
#   "<name>"           exactly that feature
# A pair takes the status of the most specific matching rule. Ties go to
# the later rule, so overrides are appended, never inserted.
frmtmb_compat_rules_tbl <- function() {
  rows <- list()
  r <- function(a, b, status, note) {
    rows[[length(rows) + 1L]] <<- data.frame(
      feature_a = a, feature_b = b, status = status, note = note,
      stringsAsFactors = FALSE)
    invisible(NULL)
  }

  ## default: nothing is assumed to work -------------------------------
  r("*", "*", "untested",
    "No rule covers this pair. The combination is not exercised by the test suite; it may work, but nothing checks it.")

  ## kind-level defaults ----------------------------------------------
  r("kind:family", "kind:covstruct", "works",
    "Covariance structures act on the linear predictor, so they are independent of the response distribution.")
  r("kind:family", "kind:special", "works",
    "Predictor specials build design columns before the family sees them.")
  r("kind:family", "kind:mode", "works",
    "Estimation modes reparameterize the outer problem and do not read the family.")
  r("kind:family", "kind:method", "conditional",
    "Depends on which post-fit ingredients the family supplies (CDF, simulator, variance function).")
  r("kind:family", "kind:aterm", "untested",
    "Addition terms are family-sensitive. See the specific rules.")
  r("kind:covstruct", "kind:covstruct", "works",
    "Separate bar terms give separate covariance blocks; the blocks do not interact.")
  r("kind:covstruct", "kind:special", "works",
    "Smooths and Gaussian processes enter as their own blocks alongside user bar terms.")
  r("kind:special", "kind:special", "works",
    "Specials are additive terms and compose freely.")
  r("kind:aterm", "kind:aterm", "works",
    "Addition terms on one response combine; each contributes its own factor to the likelihood.")
  r("kind:aterm", "kind:covstruct", "works",
    "Addition terms change the likelihood, covariance structures change the predictor.")
  r("kind:aterm", "kind:special", "works",
    "Addition terms and predictor specials are independent.")
  r("kind:mode", "kind:mode", "untested",
    "Mode pairs need a specific rule; see below.")
  r("kind:mode", "kind:covstruct", "untested",
    "See the quadrature rules; the other modes are structure-agnostic.")
  r("kind:mode", "kind:special", "untested", "See the specific rules.")
  r("kind:mode", "kind:aterm", "untested", "See the specific rules.")
  r("kind:mode", "kind:structure", "untested", "See the specific rules.")
  r("kind:structure", "kind:method", "untested", "See the specific rules.")

  ## modes that compose with everything -------------------------------
  r("verbose", "*", "works",
    "verbose only prints progress. It does not change the model or the estimates.")
  # sparse_x and autoscale are claimed only where the model surface is
  # exercised. Structures and post-fit methods keep the untested
  # default until something checks them.
  for (m in c("sparse_x", "autoscale")) {
    note <- if (m == "sparse_x") {
      "sparse_x changes how the fixed-effect design is stored. The fit is numerically the same."
    } else {
      "autoscale rescales the outer parameters for the optimizer and unscales the result."
    }
    for (k in c("family", "covstruct", "special", "aterm", "mode")) {
      r(m, paste0("kind:", k), "works", note)
    }
  }
  r("sparse_x", "mi()", "works",
    "Verified: the na.pass model frame that mi() needs and the sparse design agree. Both fits give the same estimates.")
  r("autoscale", "se()", "works",
    "Verified: known standard errors survive rescaling.")
  r("sparse_x", "s()", "works",
    "Verified: smooth basis columns join the sparse fixed-effect design.")

  ## REML --------------------------------------------------------------
  r("REML", "quadrature", "refused",
    "Refused: quadrature already marginalizes the random effects, so there is no inner problem left for REML to integrate.")
  r("REML", "profile", "refused",
    "Refused: profile = TRUE and REML both remove the fixed effects from the outer problem.")
  r("REML", "bounds", "refused",
    "Refused: under REML the fixed effects leave the outer parameter vector, so bounds naming them are rejected as unknown parameters.")
  r("REML", "priors", "conditional",
    "Priors on fixed effects are accepted under REML while bounds on the same parameters are refused. The two surfaces are inconsistent. Priors on variance parameters behave as expected.")
  r("REML", "hypothesis_profile", "refused",
    "Refused: profile likelihood tests need an ML fit.")
  r("REML", "kind:structure", "untested", "See the specific rules.")
  r("REML", "mixture", "broken",
    "Broken: mixture() under REML either stops with 'NA/NaN gradient evaluation' or returns a fit with a gradient near 1e9. There is no guard. Use REML = FALSE.")
  r("REML", "mixture_mvn", "untested",
    "Not exercised. Treat it like mixture() and prefer REML = FALSE.")
  r("REML", "mvbf", "works", "")
  r("REML", "rescor", "works", "")
  r("REML", "nl", "works", "Verified by a tiny fit.")
  r("REML", "|ID|", "works", "")
  r("REML", "mi()", "works", "Verified by a tiny fit.")
  r("REML", "cens()", "works", "Verified by a tiny fit.")
  r("REML", "trunc()", "works", "Verified by a tiny fit.")
  r("REML", "se()", "works", "Verified by a tiny fit.")
  r("REML", "gp_pred()", "works", "Verified by a tiny fit.")
  r("REML", "kind:covstruct", "works",
    "REML integrates the fixed effects and leaves the covariance blocks alone.")

  ## quadrature --------------------------------------------------------
  r("quadrature", "group:wide_blocks", "refused",
    "Refused: quadrature marginalizes one scalar random intercept at a time. Every block must be a dimension-1 us, diag, or homdiag term.")
  r("quadrature", "group:quadrature_blocks", "conditional",
    "Allowed only when the block is one-dimensional, that is a scalar random intercept. Correlated slopes are refused.")
  r("quadrature", "mi()", "refused",
    "Refused: the imputed values are themselves random effects that quadrature cannot marginalize.")
  r("quadrature", "profile", "refused",
    "Refused: profile = TRUE cannot be combined with quadrature = TRUE.")
  r("quadrature", "s()", "refused",
    "Refused in practice: a smooth is a wide random-effect block, so the scalar-intercept guard rejects it.")
  r("quadrature", "t2()", "refused",
    "Refused in practice: a smooth is a wide random-effect block.")
  r("quadrature", "gp_pred()", "refused",
    "Refused in practice: gp() builds a wide block, so the scalar-intercept guard rejects it.")
  r("quadrature", "trunc()", "broken",
    "BROKEN. The truncation normalizer is dropped from the marginalized objective. The fit returns convergence 0 with slope estimates collapsed to zero, next to a correct 0.93 from ML, REML, and profile. A large-gradient warning is the only signal. Do not combine these; use REML or plain ML for truncated responses.")
  r("quadrature", "cens()", "works",
    "Verified: agrees with the REML and ML fits of the same censored model.")
  r("quadrature", "weights()", "works", "Verified by a tiny fit.")
  r("quadrature", "se()", "works", "Verified by a tiny fit.")
  r("quadrature", "bounds", "works", "Verified by a tiny fit.")
  r("quadrature", "priors", "works", "Verified by a tiny fit.")
  r("quadrature", "mo()", "untested", "")
  r("quadrature", "nl", "works",
    "Verified by a tiny fit, with the random effect on one nonlinear parameter.")
  r("quadrature", "group:ordinal", "works", "Verified by a tiny fit.")
  r("quadrature", "mixture", "broken",
    "Broken: the fit runs and reports a gradient near 1e11. There is no guard.")
  r("quadrature", "mvbf", "untested", "")
  r("quadrature", "rescor", "untested", "")
  r("quadrature", "mixture_mvn", "untested", "")

  ## profile -----------------------------------------------------------
  r("profile", "bounds", "refused",
    "Refused: with the fixed effects profiled out they leave the outer parameter vector, so bounds naming them are rejected. This pair was once accepted and the bounds were then applied to the wrong parameters.")
  r("profile", "confint_profile", "refused",
    "Refused: profile and uniroot confidence intervals cannot address the profiled-out fixed effects. The error names the remaining parameters.")
  r("profile", "hypothesis_profile", "refused",
    "Refused: profile likelihood hypothesis tests need a fit without frmtmb_control(profile = TRUE).")
  r("profile", "priors", "conditional",
    "Priors on the profiled fixed effects are accepted, while bounds on the same parameters are refused. Treat priors under profile with care.")
  r("profile", "mixture", "broken",
    "Broken: stops with 'NA/NaN gradient evaluation' rather than a clear refusal.")
  r("profile", "mi()", "works", "Verified by a tiny fit.")
  r("profile", "cens()", "works", "Verified by a tiny fit.")
  r("profile", "trunc()", "works",
    "Verified: agrees with the ML and REML fits, unlike quadrature.")
  r("profile", "se()", "works", "Verified by a tiny fit.")
  r("profile", "weights()", "works", "Verified by a tiny fit.")
  r("profile", "kind:covstruct", "untested",
    "Profiling touches the fixed effects only, but the covariance blocks are not separately exercised.")

  ## bounds and priors -------------------------------------------------
  r("bounds", "priors", "works",
    "set_prior() carries bounds of its own; the lower and upper arguments set them directly.")
  r("bounds", "kind:aterm", "untested", "")
  r("priors", "kind:aterm", "untested", "")

  ## cens() and trunc() ------------------------------------------------
  r("cens()", "kind:family", "refused",
    "Refused: cens() needs a family with an AD log-CDF.")
  r("trunc()", "kind:family", "refused",
    "Refused: trunc() needs a family with an AD log-CDF.")
  r("trunc()", "group:cdf", "works",
    "Supported through the whole surface: the likelihood, fitted(), predict(), simulate(), and residuals().")
  r("cens()", "group:cdf_continuous", "works",
    "Supported through the whole surface; one-step-ahead residuals cover the uncensored rows only, see the residuals_osa rule.")
  r("cens()", "poisson", "refused",
    "Refused: censoring is not supported for discrete families yet, even though poisson carries a CDF.")
  r("cens()", "group:discrete", "refused",
    "Refused: censoring is not supported for discrete families yet.")
  r("trunc()", "poisson", "conditional",
    "Discrete truncation needs a lower bound of at least 1; trunc(lb = 0) is not truncation and is refused.")
  r("cens()", "trunc()", "works",
    "Both restrict the same response; the truncated likelihood is renormalized and then censored.")
  r("cens()", "weights()", "works",
    "Verified: the case weight multiplies the censored contribution.")
  r("trunc()", "weights()", "works", "")
  r("cens()", "se()", "untested", "")
  r("cens()", "mixture", "refused",
    "Refused: mixture() has no CDF, so the CDF guard rejects cens().")
  r("trunc()", "mixture", "refused", "Refused: mixture() has no CDF.")
  r("cens()", "mixture_mvn", "refused",
    "Refused: mixture_mvn() has no CDF.")
  r("trunc()", "mixture_mvn", "refused",
    "Refused: mixture_mvn() has no CDF.")
  r("cens()", "group:ordinal", "refused",
    "Refused: ordinal families carry no AD log-CDF over the response scale.")
  r("trunc()", "group:ordinal", "refused",
    "Refused: ordinal families carry no AD log-CDF over the response scale.")

  ## se() ---------------------------------------------------------------
  r("se()", "kind:family", "refused",
    "Refused: known standard errors are added to the residual variance, which only the gaussian and student families have.")
  r("se()", "group:gaussian_like", "works",
    "Supported. se(x, sigma = TRUE) keeps the estimated residual SD alongside the known one.")

  ## mi() ----------------------------------------------------------------
  r("mi()", "kind:family", "refused",
    "Refused: an imputation model must be gaussian or student.")
  r("mi()", "group:gaussian_like", "works",
    "Supported. The missing values become latent parameters of the imputation model.")
  r("mi()", "cens()", "refused",
    "Refused: mi() cannot be combined with cens(), trunc(), or se() on the same response.")
  r("mi()", "trunc()", "refused",
    "Refused: mi() cannot be combined with cens(), trunc(), or se() on the same response.")
  r("mi()", "se()", "refused",
    "Refused: mi() cannot be combined with cens(), trunc(), or se() on the same response.")
  r("mi()", "rescor", "refused",
    "Refused: mi() cannot be combined with rescor = TRUE.")
  r("mi()", "mvbf", "works",
    "Required, in fact: mi() needs a second bf() supplying the imputation model for the variable.")
  r("mi()", "mi_pred()", "conditional",
    "The predictor mi(x) needs a matching response model x | mi() ~ ... in the same mvbf(). A variable may not impute itself.")
  r("mi()", "weights()", "untested", "")

  ## weights() and trials() ---------------------------------------------
  r("weights()", "kind:family", "works",
    "Case weights multiply each observation's log-likelihood contribution, whatever the family.")
  r("trials()", "kind:family", "untested",
    "trials() is meaningful only for the binomial-type families.")
  r("trials()", "group:trials_families", "works", "")
  r("trials()", "multinomial", "conditional",
    "Required. The row sums of the response matrix must equal the trials.")

  ## vint() and vreal() --------------------------------------------------
  r("vint()", "*", "untested",
    "vint() passes integer covariates to a custom family; only the custom-family path uses them.")
  r("vreal()", "*", "untested",
    "vreal() passes real covariates to a custom family; only the custom-family path uses them.")
  r("vint()", "vreal()", "works", "Both may appear on one response.")
  r("vint()", "mvbf", "untested",
    "Custom-family covariates inside a multivariate model are not exercised.")
  r("vreal()", "mvbf", "untested",
    "Custom-family covariates inside a multivariate model are not exercised.")

  ## rescor ---------------------------------------------------------------
  r("rescor", "kind:family", "refused",
    "Refused: rescor = TRUE requires every response to be gaussian.")
  r("rescor", "gaussian", "works", "")
  r("rescor", "cens()", "refused",
    "Refused. This pair was once accepted with the censoring silently dropped.")
  r("rescor", "trunc()", "refused",
    "Refused. This pair was once accepted with the truncation silently dropped.")
  r("rescor", "se()", "refused", "Refused.")
  r("rescor", "weights()", "refused", "Refused.")
  r("rescor", "mvbf", "works",
    "rescor = TRUE is only meaningful inside mvbf(); set_rescor() on a single formula is refused.")
  r("rescor", "kind:method", "untested", "")
  r("rescor", "fitted", "works", "")
  r("rescor", "predict", "works", "")
  r("rescor", "simulate", "refused",
    "Refused: simulate() is not supported for multivariate fits yet.")
  r("rescor", "residuals_osa", "refused",
    "Refused: residuals() is not supported for multivariate fits yet.")
  r("rescor", "emmeans", "refused",
    "Refused: emmeans support is univariate-only for now.")

  ## mvbf ------------------------------------------------------------------
  r("mvbf", "kind:method", "refused",
    "Refused: the post-fit methods below are univariate-only for now.")
  r("mvbf", "fitted", "works", "")
  r("mvbf", "predict", "works", "")
  r("mvbf", "kind:family", "works",
    "Each response carries its own family unless rescor = TRUE.")
  r("mvbf", "kind:aterm", "works",
    "Addition terms are per response.")
  r("mvbf", "|ID|", "works",
    "|ID| is the reason mvbf() exists in this grammar: it keys correlated random effects across responses.")

  ## |ID| -------------------------------------------------------------------
  r("|ID|", "kind:covstruct", "refused",
    "Refused: |ID| correlation is only supported for default (us) random-effect terms.")
  r("|ID|", "us", "works", "")
  r("|ID|", "rr", "refused",
    "Refused: rr() terms cannot share an |ID| key.")

  ## nl ----------------------------------------------------------------------
  r("nl", "kind:family", "refused",
    "Refused: nl = TRUE requires a family with a single mu location parameter.")
  r("nl", "gaussian", "works", "")
  r("nl", "kind:method", "untested", "")
  r("nl", "fitted", "works", "")
  r("nl", "predict", "conditional",
    "Point predictions work. se.fit is not supported for the nonlinear predictor; request a nonlinear parameter with dpar instead.")
  r("nl", "cens()", "works", "Verified by a tiny fit.")
  r("nl", "kind:covstruct", "works",
    "A nonlinear parameter may carry its own random effects.")

  ## mixture -----------------------------------------------------------------
  r("mixture", "mvbf", "untested", "")
  r("mixture", "rescor", "refused",
    "Refused: rescor = TRUE requires gaussian responses, and a mixture is not one.")
  r("mixture", "mi()", "refused",
    "Refused: mi() on the mixture response is not supported.")
  r("mixture", "se()", "refused",
    "Refused: se() is supported for gaussian and student families only.")
  r("mixture", "weights()", "works", "Verified by a tiny fit.")
  r("mixture", "simulate", "conditional",
    "Works only when every component family has a simulator.")
  r("mixture", "frm_sample", "conditional",
    "Mixture posteriors are multimodal. Sample with init = \"random\" rather than the mode-anchored default.")
  r("mixture", "kind:covstruct", "untested", "")
  r("mixture_mvn", "simulate", "refused",
    "Refused: mixture_mvn() has no simulator yet.")
  r("mixture_mvn", "mvbf", "refused",
    "Refused: the family already takes a matrix response.")
  r("mixture_mvn", "kind:aterm", "untested", "")

  ## specials ------------------------------------------------------------------
  r("mo()", "kind:family", "conditional",
    "The monotonic variable must be an ordered factor or an integer with at least 3 categories.")
  r("mo()", "mo()", "refused",
    "Refused: mo() cannot interact with another mo() or mi() term.")
  r("mo()", "mi_pred()", "refused",
    "Refused: mo() cannot interact with another mo() or mi() term.")
  r("mo()", "kind:mode", "untested", "")
  r("mo()", "REML", "works", "Verified by a tiny fit.")
  r("mo()", "profile", "works", "Verified by a tiny fit.")
  r("mo()", "predict", "conditional",
    "New data must stay inside the fitted category range; unknown categories are refused.")
  r("cs_pred()", "kind:family", "refused",
    "Refused: cs() needs an sratio, cratio, or acat family.")
  r("cs_pred()", "group:ordinal_cs", "works", "")
  r("cs_pred()", "cumulative", "refused",
    "Refused: category-specific effects are not identified under the cumulative parameterization.")
  r("gp_pred()", "kind:mode", "untested", "")
  r("gp_pred()", "kind:family", "works", "")
  r("s()", "kind:mode", "untested", "")
  r("s()", "kind:family", "works", "")
  r("t2()", "kind:family", "works", "")

  ## covariance structures -------------------------------------------------------
  r("group:positional", "kind:special", "conditional",
    "Positional structures take the lag from the level order, not from the level value.")
  r("ar1", "*", "conditional",
    "ar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns. Use ou() over num_factor() for irregular spacing.")
  r("hetar1", "*", "conditional",
    "hetar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns.")
  r("ou", "*", "conditional",
    "ou() is the irregular-spacing structure. Build the levels with num_factor() so the metric distance is recoverable.")
  r("group:spatial", "*", "conditional",
    "Spatial structures need coordinates built with num_factor(x, y).")
  r("gr_prec", "*", "conditional",
    "gr(prec = Q) supports intercept-only terms: (1 | gr(g, prec = Q)).")
  r("gr_cov", "*", "conditional",
    "gr(cov = A) accepts correlated slopes; the block covariance is the Kronecker product of A and the term covariance. A needs dimnames covering every grouping level.")
  r("equalto", "*", "conditional",
    "equalto(x + 0 | g, V) fixes the term covariance to V, which must be square and match the term dimension.")
  r("rr", "*", "conditional",
    "rr() gives a reduced-rank block; the rank d must not exceed the term dimension.")
  r("rr", "sparse_x", "works",
    "Verified: identical estimates with and without a sparse fixed-effect design.")
  r("rr", "REML", "works", "Verified by a tiny fit.")
  r("smooth", "*", "conditional",
    "smooth is the internal structure behind s() and t2(); it is not written directly in a formula.")
  r("group:latent_gp", "*", "conditional",
    "gp() and hsgp() are predictor specials, not bar terms. Write gp(x), not (gp(x) | g).")
  r("gp_pred()", "gp_pred()", "conditional",
    "gp() takes 1 to 3 variables. The arguments k, c, and iso are evaluated in the formula environment, and c may be a vector with one entry per dimension.")

  ## post-fit methods --------------------------------------------------------------
  r("simulate", "kind:family", "works",
    "The family supplies a simulator.")
  r("simulate", "group:no_simulator", "refused",
    "Refused: this family has no simulator yet.")
  r("residuals_osa", "kind:family", "conditional",
    "One-step-ahead residuals need the family to register its observation through OBS().")
  r("residuals_osa", "group:ordinal", "works",
    "oneStepGeneric over the discrete support 1..K; the result is a randomized quantile residual and matches the analytic one to 1e-13.")
  r("residuals_osa", "cens()", "conditional",
    "Censored rows return NA: what is observed there is an event, not a value, so it carries no one-step CDF. The uncensored rows get residuals conditional on the censoring events, which needs one censoring point per side (type-I censoring). Row-varying censoring points and interval censoring are refused.")
  r("residuals_osa", "trunc()", "conditional",
    "Supported for constant truncation bounds. Row-varying bounds are refused, because the one-step-ahead transform is not defined across changing support.")
  r("residuals_osa", "weights()", "conditional",
    "Runs, but the residuals ignore the case weights. Treat them as unweighted.")
  r("residuals_osa", "kind:structure", "untested", "")
  r("emmeans", "kind:family", "conditional",
    "Univariate fits only, and the mu predictor must be linear.")
  r("emmeans", "nl", "refused",
    "Refused: emmeans support needs a linear mu predictor.")
  r("frm_sample", "*", "conditional",
    "Chains start jittered around the fitted mode. Use init_jitter to widen the spread, or init = \"random\" for a multimodal posterior.")
  r("confint_profile", "kind:mode", "untested", "")
  r("hypothesis_profile", "kind:mode", "untested", "")
  r("predict", "kind:family", "conditional",
    "Rank-deficient designs drop aliased columns at fit time. New data that is not estimable from the retained columns predicts NA and warns.")
  r("fitted", "kind:family", "conditional",
    "Needs a family with a mean function.")

  ## formula grammar --------------------------------------------------------------
  r("bar_crossing", "*", "refused",
    "Refused: a bar term crossed with * or : (as in x * (1 | g)) is not a random-effect specification (lme4#196). Write the crossing inside the bar: (x | g). This spelling was once accepted with the crossing silently dropped.")
  r("call_group", "*", "works",
    "Call-valued grouping factors are supported: (1 | factor(x)) and (1 | interaction(a, b)) both build the grouping factor from the model frame.")
  r("double_bar", "*", "works",
    "(x || g) gives uncorrelated terms. With a factor on the left, (f || g) routes to diag, that is one independent effect per factor level.")
  r("double_bar", "diag", "works",
    "diag is the structure (f || g) resolves to.")
  r("bar_crossing", "kind:covstruct", "refused",
    "Refused for every covariance structure. The crossing must go inside the bar.")

  ## the family does not change a covariance structure's own condition ---------
  for (cs in c("ar1", "hetar1", "ou", "exp", "gau", "mat", "rr",
               "equalto", "gr_cov", "gr_prec", "smooth", "gp",
               "hsgp")) {
    r(cs, "kind:family", "works",
      "The covariance structure acts on the linear predictor, so the response distribution does not change its conditions.")
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------- resolution

frmtmb_compat_statuses <- c("works", "conditional", "refused",
                            "broken", "untested")

# Specificity of one pattern. Higher beats lower.
compat_spec <- function(pat) {
  ifelse(pat == "*", 0L,
         ifelse(startsWith(pat, "kind:"), 1L,
                ifelse(startsWith(pat, "group:"), 2L, 3L)))
}

# Does `pat` cover the feature at position i of the feature table?
compat_match <- function(pat, name, kind) {
  if (pat == "*") return(rep(TRUE, length(name)))
  if (startsWith(pat, "kind:")) return(kind == substring(pat, 6L))
  if (startsWith(pat, "group:")) {
    g <- frmtmb_compat_groups_lst[[substring(pat, 7L)]]
    return(name %in% g)
  }
  name == pat
}

#' Feature metadata for the compatibility registry
#'
#' The vocabulary the compatibility registry talks about: every family,
#' addition term, covariance structure, predictor special, estimation
#' mode, model structure, and post-fit method that has a declared
#' compatibility status.
#'
#' `name` is how the feature is written in a formula or a call. `key`
#' is the identifier the package uses internally, which is what lets
#' the tests check the registry against [frm()]'s real vocabulary.
#' Three specials share a name with a covariance structure, so they
#' carry the display names `mi_pred()`, `gp_pred()`, and `cs_pred()`.
#'
#' @return A data frame with columns `name`, `key`, and `kind`.
#' @seealso [frm_compat()], [frm_compat_rules()]
#' @examples
#' head(frm_compat_features())
#' table(frm_compat_features()$kind)
#' @export
frm_compat_features <- function() {
  frmtmb_compat_features_tbl()
}

#' Compatibility rules, before resolution
#'
#' The registry stores rules, not one row per pair. A rule side is one
#' of four patterns, from least to most specific: `"*"` for any
#' feature, `"kind:<kind>"` for a whole kind, `"group:<group>"` for a
#' named feature set, or a bare feature name. A pair takes the status
#' of the most specific matching rule; ties go to the later rule.
#'
#' Use [frm_compat()] to read the resolved answer for a pair. Use this
#' function to see which rule is doing the work.
#'
#' @return A data frame with columns `feature_a`, `feature_b`,
#'   `status`, and `note`.
#' @seealso [frm_compat()], [frm_compat_features()]
#' @examples
#' rules <- frm_compat_rules()
#' nrow(rules)
#' subset(rules, status == "broken")[, c("feature_a", "feature_b")]
#' @export
frm_compat_rules <- function() {
  frmtmb_compat_rules_tbl()
}

#' Query the feature compatibility registry
#'
#' Gives the declared status of a pair of package features. The status
#' is one of:
#'
#' \describe{
#'   \item{`works`}{The combination is supported and exercised.}
#'   \item{`conditional`}{Supported, but the note states a condition
#'     that the combination must meet.}
#'   \item{`refused`}{[frm()] or the post-fit method stops with an
#'     error. The refusal is deliberate.}
#'   \item{`broken`}{The combination is accepted but the result is
#'     wrong, or it fails with an error that does not explain itself.
#'     Avoid the pair. The note gives the evidence.}
#'   \item{`untested`}{Nothing checks this pair. It may work. Treat a
#'     silent success as unverified, not as support.}
#' }
#'
#' The last status is the point of the registry. A guard that does not
#' exist looks exactly like a guard that passed, so the absence of an
#' error was never evidence of support.
#'
#' @param feature_a,feature_b Feature names, as given by
#'   [frm_compat_features()]. Supply both for one pair, one for every
#'   pair involving that feature, or neither for the whole table.
#' @param status Optional character vector; keep only these statuses.
#' @return A data frame with columns `feature_a`, `kind_a`,
#'   `feature_b`, `kind_b`, `status`, and `note`. Family pairs are
#'   omitted, because a model carries one family.
#' @seealso [frm_compat_rules()], [frm_compat_features()]
#' @examples
#' # one pair
#' frm_compat("rescor", "cens()")
#'
#' # everything known about truncation
#' frm_compat("trunc()", status = c("refused", "broken"))
#'
#' # the pairs to avoid
#' frm_compat(status = "broken")[, 1:5]
#' @export
frm_compat <- function(feature_a = NULL, feature_b = NULL,
                       status = NULL) {
  ft <- frmtmb_compat_features_tbl()
  for (nm in c(feature_a, feature_b)) {
    if (!nm %in% ft$name) {
      stop("Unknown feature: '", nm, "'. See frm_compat_features().",
           call. = FALSE)
    }
  }
  if (!is.null(status)) {
    bad <- setdiff(status, frmtmb_compat_statuses)
    if (length(bad)) {
      stop("Unknown status: ", paste(bad, collapse = ", "),
           ". Statuses are: ",
           paste(frmtmb_compat_statuses, collapse = ", "),
           call. = FALSE)
    }
  }

  # unordered pairs of distinct features, minus family x family: a
  # model carries exactly one family, so those pairs mean nothing
  n <- nrow(ft)
  ij <- utils::combn(n, 2L)
  keep <- !(ft$kind[ij[1, ]] == "family" & ft$kind[ij[2, ]] == "family")
  ij <- ij[, keep, drop = FALSE]
  pairs <- data.frame(
    feature_a = ft$name[ij[1, ]], kind_a = ft$kind[ij[1, ]],
    feature_b = ft$name[ij[2, ]], kind_b = ft$kind[ij[2, ]],
    status = NA_character_, note = NA_character_,
    stringsAsFactors = FALSE)

  # Resolve before subsetting so a rule's specificity is judged against
  # the whole table, never against the slice the caller asked for.
  best <- rep(-1L, nrow(pairs))
  rules <- frmtmb_compat_rules_tbl()
  for (k in seq_len(nrow(rules))) {
    pa <- rules$feature_a[k]
    pb <- rules$feature_b[k]
    spec <- compat_spec(pa) + compat_spec(pb)
    # unordered: the rule may match the pair either way round
    hit <- (compat_match(pa, pairs$feature_a, pairs$kind_a) &
              compat_match(pb, pairs$feature_b, pairs$kind_b)) |
      (compat_match(pa, pairs$feature_b, pairs$kind_b) &
         compat_match(pb, pairs$feature_a, pairs$kind_a))
    # >= so a later rule of equal specificity overrides an earlier one
    hit <- hit & spec >= best
    pairs$status[hit] <- rules$status[k]
    pairs$note[hit] <- rules$note[k]
    best[hit] <- spec
  }

  sel <- rep(TRUE, nrow(pairs))
  if (!is.null(feature_a) && !is.null(feature_b)) {
    sel <- (pairs$feature_a == feature_a & pairs$feature_b == feature_b) |
      (pairs$feature_a == feature_b & pairs$feature_b == feature_a)
  } else {
    one <- feature_a %||% feature_b
    if (!is.null(one)) {
      sel <- pairs$feature_a == one | pairs$feature_b == one
    }
  }
  if (!is.null(status)) sel <- sel & pairs$status %in% status
  out <- pairs[sel, , drop = FALSE]

  # read the queried feature first, so a one-feature slice reads down a
  # single column
  one <- if (is.null(feature_b)) feature_a else feature_b
  if (!is.null(one) && (is.null(feature_a) || is.null(feature_b))) {
    flip <- out$feature_b == one
    if (any(flip)) {
      out[flip, c("feature_a", "kind_a", "feature_b", "kind_b")] <-
        out[flip, c("feature_b", "kind_b", "feature_a", "kind_a")]
    }
  }
  rownames(out) <- NULL
  out
}
