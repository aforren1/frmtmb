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

# ------------------------------------------------------ contributor seam
#
# A feature that does not live in this file contributes its own rows.
# The structured families went this way while they were still here, and
# it is what let them leave with their rules attached: the matrix covers
# a family in another package without R/compat.R knowing it exists.
#
# Contributions are APPENDED, and that is the only safe direction. Rules
# of equal specificity are resolved by position, later wins (see
# frmtmb_compat_rules_tbl()), so a contributed rule may override a core
# default and a core default can never silently override a contributed
# one.
#
# The container is created here and filled by the contributors. Every
# contributor is now an outside package, which registers in its own
# .onLoad(): by then every namespace is sealed, so the collation-order
# question a top-level call inside this package would raise does not
# arise.

frmtmb_compat_contrib <- new.env(parent = emptyenv())
frmtmb_compat_contrib$features <- list()
frmtmb_compat_contrib$rules <- list()

#' Contribute to the compatibility matrix from another package
#'
#' [frm_compat()] answers what one feature does in the presence of
#' another, and it answers in three states: the pair works, the pair is
#' refused, or the pair is untested. A feature that lives outside
#' frmtmb has to be able to say the same things about itself, or the
#' matrix reports on part of the package and reads as if it reported on
#' all of it.
#'
#' `frmtmb_register_compat()` is that seam. Call it from the
#' contributing package's `.onLoad()`: by then every namespace is
#' sealed, so the collation-order question a top-level call would raise
#' does not arise. `compat_rule_builder()` is the accumulator the core's
#' own rules are written with, exported so that a contributed rule reads
#' the same as a core one.
#'
#' Contributions are APPENDED, which is the only safe direction. Rules
#' of equal specificity resolve later-wins, so a contributed rule may
#' override a core default and a core default can never silently
#' override a contributed one.
#'
#' @section Status vocabulary:
#' Every rule declares one of three states, and the third is the reason
#' the registry exists: a guard that does not exist looks exactly like a
#' guard that passed, so "no error" was never evidence of support.
#' \describe{
#'   \item{`"works"`}{The pair is supported and exercised by a test.}
#'   \item{`"refused"`}{The pair is rejected with a message. Say in the
#'     note what the user should do instead.}
#'   \item{`"untested"`}{Nothing is known. Neither a promise nor a
#'     refusal.}
#' }
#'
#' @param features Named character vector mapping a feature's DISPLAY
#'   name to its kind: `c("hmm" = "structure")`, `c("dec()" = "aterm")`.
#'   The display name is how the feature is written in a formula or a
#'   call, parentheses included for a callable one; the kind groups it
#'   in the printed matrix (`"family"`, `"covstruct"`, `"aterm"`,
#'   `"special"`, `"autocor"`, `"mode"`, `"structure"`).
#' @param rules A FUNCTION of no arguments returning a rule data frame,
#'   not the data frame itself, so that contributed rules are built on
#'   demand exactly as the core ones are. Build the frame with
#'   `compat_rule_builder()`: it returns a list of two functions, `r(a,
#'   b, status, note, override = FALSE)` to record one rule and
#'   `rules()` to return the accumulated frame. `a` and `b` are display
#'   names, `"*"` matches every feature, and `override = TRUE` says the
#'   rule beats a more specific one rather than losing to it.
#' @return `NULL`, invisibly. Called for the registration.
#'   `compat_rule_builder()` returns a list with elements `r` and
#'   `rules`.
#' @seealso [frm_compat()] for the matrix these fill,
#'   [frmtmb_family()] and [frmtmb_structure()] for the family-side
#'   seams a contributor usually registers alongside these, and
#'   `vignette("compatibility")`
#' @examples
#' # what a contributing package's .onLoad() does
#' contribute <- function() {
#'   b <- compat_rule_builder()
#'   b$r("wiener", "cens()", "refused",
#'       "The family supplies no lcdf, so there is no CDF to censor with.")
#'   b$r("wiener", "*", "untested",
#'       "Not exercised outside this package's own suite.")
#'   frmtmb_register_compat(features = c("wiener" = "family"),
#'                          rules = b$rules)
#' }
#' # the accumulator on its own, which is all a rule set is
#' b <- compat_rule_builder()
#' b$r("dec()", "trials()", "refused", "Different response shapes.")
#' b$rules()
#' @export
frmtmb_register_compat <- function(features = NULL, rules = NULL) {
  if (length(features)) {
    frmtmb_compat_contrib$features <-
      c(frmtmb_compat_contrib$features, list(features))
  }
  if (!is.null(rules)) {
    frmtmb_compat_contrib$rules <-
      c(frmtmb_compat_contrib$rules, list(rules))
  }
  invisible(NULL)
}

#' @rdname frmtmb_register_compat
#' @export
compat_rule_builder <- function() {
  rows <- list()
  list(
    r = function(a, b, status, note, override = FALSE) {
      rows[[length(rows) + 1L]] <<- data.frame(
        feature_a = a, feature_b = b, status = status, note = note,
        override = override, stringsAsFactors = FALSE)
      invisible(NULL)
    },
    rules = function() {
      out <- do.call(rbind, rows)
      rownames(out) <- NULL
      out
    }
  )
}

# ---------------------------------------------------------------- features

#' The feature vocabulary of the registry, one row per feature.
#'
#' `name`: how the feature is written in a formula or a call argument.
#' `key`:  the identifier the package itself uses, so tests can look the
#'         feature up in `family_registry`, `covstruct_registry`, and the
#'         parser vocabularies. Display names carry `"()"` for callable
#'         features; the key never does.
#'
#' @noRd
frmtmb_compat_features_tbl <- function() {
  contrib <- unlist(frmtmb_compat_contrib$features)
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
            "acat", "categorical", "von_mises", "cox")
  covs <- c("us", "diag", "homdiag", "cs", "ar1", "hetar1", "ou",
            "toep", "homtoep", "homcs", "exp", "gau", "mat", "rr",
            "equalto", "gr_cov", "gr_prec", "smooth", "gp", "hsgp",
            "car", "spde", "us_t", "diag_t")
  do.call(rbind, c(
    lapply(fams, f, kind = "family"),
    lapply(covs, f, kind = "covstruct"),
    lapply(c("weights()", "trials()", "cens()", "trunc()", "se()",
             "mi()", "vint()", "vreal()"), f, kind = "aterm"),
    lapply(c("s()", "t2()", "mo()", "mi_pred()", "gp_pred()",
             "cs_pred()"), f, kind = "special"),
    # R-side (within-group residual) correlation terms. They carry no
    # random effect, so they are not a covariance structure, and they
    # contribute no design column, so they are not a predictor special:
    # they replace the response's density. Hence a kind of their own.
    lapply(c("ar()", "ma()", "arma()", "cosy()", "unstr()"), f,
           kind = "autocor"),
    lapply(c("REML", "quadrature", "profile", "autoscale", "sparse_x",
             "prior", "bounds", "verbose"), f, kind = "mode"),
    lapply(c("mvbf", "rescor", "|ID|", "nl", "mixture",
             "mixture_mvn"), f, kind = "structure"),
    lapply(c("fitted", "predict", "simulate", "residuals",
             "residuals_osa", "emmeans",
             "confint_profile", "hypothesis_profile"), f,
           kind = "method"),
    # formula-grammar spellings, which have their own restrictions and
    # belong in the table even though they name no package object
    lapply(c("bar_crossing", "call_group", "double_bar", "mm()",
             "mmc()"), f, kind = "grammar"),
    # contributed last, so that the vocabulary a contributor adds cannot
    # displace a core feature's position in the pair table
    lapply(seq_along(contrib), function(i) {
      f(names(contrib)[[i]], unname(contrib[[i]]))
    })
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
          "inverse.gaussian", "cox"),
  # cens() additionally refuses discrete responses, so poisson drops out
  cdf_continuous = c("gaussian", "lognormal", "exponential", "weibull",
                     "inverse.gaussian", "cox"),
  # families whose modelled response is a distribution over categories
  # rather than a number, so fitted() returns an n x K matrix
  categorical_probs = c("cumulative", "sratio", "cratio", "acat",
                        "categorical"),
  gaussian_like = c("gaussian", "student"),
  ordinal = c("cumulative", "sratio", "cratio", "acat"),
  # cs() category-specific effects are undefined under the cumulative
  # parameterization; the sequential and adjacent-category ones take them
  ordinal_cs = c("sratio", "cratio", "acat"),
  discrete = c("poisson", "negbinomial", "nbinom1", "geometric",
               "compois", "binomial", "bernoulli", "beta_binomial",
               "zero_inflated_poisson", "zero_inflated_negbinomial",
               "zero_inflated_binomial", "hurdle_poisson"),
  no_simulator = c("tweedie", "compois", "hurdle_poisson", "cox"),
  matrix_response = c("multinomial"),
  trials_families = c("binomial", "beta_binomial",
                      "zero_inflated_binomial", "multinomial"),
  # quadrature marginalizes one scalar random effect at a time
  quadrature_blocks = c("us", "diag", "homdiag"),
  # every covariance structure quadrature cannot marginalize
  wide_blocks = c("cs", "ar1", "hetar1", "ou", "toep", "homtoep",
                  "homcs", "exp", "gau", "mat", "rr", "equalto",
                  "gr_cov", "gr_prec", "smooth", "gp", "hsgp",
                  "car", "spde"),
  # sparse Gaussian Markov random fields written as predictor specials
  gmrf = c("car", "spde"),
  # the R-side residual correlation terms, which share every guard
  autocor = c("ar()", "ma()", "arma()", "cosy()", "unstr()"),
  # the two of them brms also treats as true residual covariance (its
  # "natural residuals" families)
  autocor_families = c("gaussian", "student"),
  # structures over a metric coordinate rather than positional levels
  spatial = c("exp", "gau", "mat"),
  # positional structures: the level order sets the lag, not the value
  positional = c("ar1", "hetar1", "toep", "homtoep"),
  latent_gp = c("gp", "hsgp"),
  known_cov = c("gr_cov", "gr_prec", "equalto"),
  # the Student-t latents, gr(dist = "student"). They are the us/diag
  # blocks with a t density over the levels, so they share the us/diag
  # guards and add the ones the per-level mixing variable imposes
  student_blocks = c("us_t", "diag_t"),
  post_fit = c("fitted", "predict", "simulate", "residuals",
               "residuals_osa", "emmeans")
)

# ------------------------------------------------------------------- rules

#' The declared compatibility rules, one row per rule.
#'
#' Patterns, from least to most specific, with the specificity each one
#' scores:
#'   `"*"`                any feature                            0
#'   `"kind:<kind>"`      any feature of that kind               1
#'   `"group:<group>"`    any feature in that named group        2
#'   `"<name>"`           exactly that feature                   3
#'
#' PRECEDENCE. A rule's signature is its two sides' specificities sorted
#' descending, and two rules are compared lexicographically on that
#' signature:
#'
#'   `(3,3) > (3,2) > (3,1) > (3,0) > (2,2) > (2,1) > (2,0) > (1,1) > ...`
#'
#' Adding the two sides instead would let unrelated pattern shapes tie -
#' `"<name>" x "*"`, `"group:" x "kind:"` and `"kind:" x "group:"` all
#' sum to 3 - and file order, not specificity, would then decide. The
#' sorted pair keeps the property that matters: a rule that is strictly
#' more specific on one side and no less specific on the other always
#' wins, because sorting preserves that domination.
#'
#' Two rules can therefore tie only when their signatures are identical.
#' The later rule wins those, so an override is appended, never
#' inserted. A tie between rules that DISAGREE about the status is a
#' registry defect, not an override: `frmtmb_compat_validate()`
#' (asserted by test-compat.R) errors on it unless the winning rule is
#' marked `override = TRUE` and says in its note what it is overriding.
#'
#' @noRd
frmtmb_compat_rules_tbl <- function() {
  rows <- list()
  r <- function(a, b, status, note, override = FALSE) {
    rows[[length(rows) + 1L]] <<- data.frame(
      feature_a = a, feature_b = b, status = status, note = note,
      override = override, stringsAsFactors = FALSE)
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
  # No kind-level default for family x aterm: addition terms are
  # family-sensitive, and every declared aterm states its own family
  # rule below, so such a default could never win a pair. A new aterm
  # added without one falls to the untested default, which is what the
  # rule used to say anyway.
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
  # No kind-level default for mode x mode either: every pair of the
  # eight declared modes is already named below, either by an explicit
  # rule or by the verbose / sparse_x / autoscale sweeps.
  r("kind:mode", "kind:covstruct", "untested",
    "See the quadrature rules; the other modes are structure-agnostic.")
  r("kind:mode", "kind:special", "untested", "See the specific rules.")
  r("kind:mode", "kind:aterm", "untested", "See the specific rules.")
  r("kind:mode", "kind:structure", "untested", "See the specific rules.")
  r("kind:structure", "kind:method", "untested", "See the specific rules.")

  ## blanket feature defaults, weakest first --------------------------
  #
  # A "<name> x *" rule is one feature's default against everything not
  # otherwise named. Two of them meeting have the same signature, so
  # only file order separates them; the order below is the authority
  # order, and every rule that relies on it carries override = TRUE:
  #
  #   permissive defaults (here) < untested defaults (vint, vreal)
  #     < a structure's own condition < an outright refusal
  #
  # so a refusal is never lost to a permissive default, and a condition
  # is never lost to "works".
  r("verbose", "*", "works",
    "verbose only prints progress. It does not change the model or the estimates.")
  r("call_group", "*", "works",
    "Call-valued grouping factors are supported: (1 | factor(x)) and (1 | interaction(a, b)) both build the grouping factor from the model frame.")
  r("double_bar", "*", "works",
    "(x || g) gives uncorrelated terms. With a factor on the left, (f || g) routes to diag, that is one independent effect per factor level.")
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
  r("REML", "prior", "conditional",
    "Priors on fixed effects are accepted under REML while bounds on the same parameters are refused. The two surfaces are inconsistent. Priors on variance parameters behave as expected.")
  r("REML", "hypothesis_profile", "refused",
    "Refused: profile likelihood tests need an ML fit.")
  # Every model structure carries its own REML rule below, so a
  # kind-level REML x structure default could never win a pair.
  r("REML", "mixture", "refused",
    "Refused: a mixture likelihood is invariant to permuting its components, so it is multimodal in the fixed effects REML integrates out. The inner Laplace solve has no single mode to expand about; the fit used to stop with 'NA/NaN gradient evaluation' or report a gradient near 1e9.")
  r("REML", "mixture_mvn", "refused",
    "Refused for the same reason as mixture(): the classes are exchangeable, so the restricted likelihood is not defined.")
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
    "Allowed only when the block is one-dimensional, that is a scalar random intercept. Correlated slopes are refused. Several such blocks are fine, nested ones included: (1 | ga/gb) becomes an iterated one-dimensional integral.")
  r("quadrature", "group:post_fit", "works",
    "The marginalized objective carries no conditional modes, so they are recovered from the inner solve of the Laplace objective at the quadrature optimum. Verified: they match glmer(nAGQ = 25)'s modes to 1e-4, and ranef(), fitted(), predict() and residuals() are NA-free.")
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
  r("quadrature", "trunc()", "refused",
    "Refused: the truncation normalizer is log(F(ub) - F(lb)) over plain CDFs, and the Gauss-Kronrod nodes reach random-effect values where that difference underflows to exactly zero while the density is still representable. The integrand is then +Inf and the marginalized objective is -Inf, at the Laplace optimum as well as at the starting values, so the fit used to report logLik = +Inf as converged. Laplace stays near the mode and is unaffected: use quadrature = FALSE, REML, or profile for truncated responses.")
  r("quadrature", "cens()", "works",
    "Verified: agrees with the REML and ML fits of the same censored model.")
  r("quadrature", "weights()", "works", "Verified by a tiny fit.")
  r("quadrature", "se()", "works", "Verified by a tiny fit.")
  r("quadrature", "bounds", "works", "Verified by a tiny fit.")
  r("quadrature", "prior", "works", "Verified by a tiny fit.")
  r("quadrature", "mo()", "untested", "")
  r("quadrature", "nl", "works",
    "Verified by a tiny fit, with the random effect on one nonlinear parameter.")
  r("quadrature", "group:ordinal", "works", "Verified by a tiny fit.")
  r("quadrature", "mixture", "conditional",
    "Exact when the per-group integrand is univariate (one scalar random intercept, in one class), which is the case the exactness test pins down. A mixture whose fit collapses a mixing weight to about exp(-35) degenerates the rescaled integrand and leaves a gradient near 1e14; the fit is not silent, the large-gradient and non-convergence warnings both fire. Check the gradient before trusting a mixture quadrature fit.")
  r("quadrature", "mvbf", "untested", "")
  r("quadrature", "rescor", "untested", "")
  r("quadrature", "mixture_mvn", "refused",
    "Refused in practice: mixture_mvn() clusters a matrix response and carries no scalar random-intercept block, so the scalar-intercept guard rejects it.")

  ## profile -----------------------------------------------------------
  r("profile", "bounds", "refused",
    "Refused: with the fixed effects profiled out they leave the outer parameter vector, so bounds naming them are rejected. This pair was once accepted and the bounds were then applied to the wrong parameters.")
  r("profile", "confint_profile", "refused",
    "Refused: profile and uniroot confidence intervals cannot address the profiled-out fixed effects. The error names the remaining parameters.")
  r("profile", "hypothesis_profile", "refused",
    "Refused: profile likelihood hypothesis tests need a fit without frmtmb_control(profile = TRUE).")
  r("profile", "prior", "conditional",
    "Priors on the profiled fixed effects are accepted, while bounds on the same parameters are refused. Treat priors under profile with care.")
  r("profile", "mixture", "refused",
    "Refused: profiling moves the fixed effects into the inner Laplace problem, and a mixture likelihood is multimodal in them. The fit used to stop with 'NA/NaN gradient evaluation' rather than a clear refusal.")
  r("profile", "mixture_mvn", "refused",
    "Refused for the same reason as mixture(): the inner problem would be multimodal.")
  r("profile", "mi()", "works", "Verified by a tiny fit.")
  r("profile", "cens()", "works", "Verified by a tiny fit.")
  r("profile", "trunc()", "works",
    "Verified: agrees with the ML and REML fits, unlike quadrature.")
  r("profile", "se()", "works", "Verified by a tiny fit.")
  r("profile", "weights()", "works", "Verified by a tiny fit.")
  r("profile", "kind:covstruct", "untested",
    "Profiling touches the fixed effects only, but the covariance blocks are not separately exercised.")

  ## bounds and priors -------------------------------------------------
  r("bounds", "prior", "works",
    "set_prior() carries bounds of its own; the lower and upper arguments set them directly.")
  r("bounds", "kind:aterm", "untested", "")
  r("prior", "kind:aterm", "untested", "")

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
    "The response is truncated FIRST and censored inside the window: a right-censored row contributes (F(ub) - F(y)) / Z and a left-censored one (F(y) - F(lb)) / Z, with Z the window mass (F(lb - 1) for a discrete lower bound). Verified against a hand-rolled likelihood for gaussian and for the discrete composed form, and it is the likelihood simulate(censored = TRUE) draws from. Through v0.25 only the normalizer was windowed, which censored the UNtruncated variable and biased the residual sd upward.")
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
  r("cens()", "cox", "works",
    "The reason the family exists. An event contributes the density log h0(t) + eta - H0(t) exp(eta) and a censored row the survivor function, which is what cox()'s log-density and log-CDF already are, so right, left and interval censoring all run through the ordinary cens() machinery. Verified against a hand-rolled M-spline PH likelihood exactly and against survival::coxph() to 2e-2 on the log hazard ratio.")
  r("cens()", "categorical", "refused",
    "Refused: a nominal response carries no order, so it has no CDF for a censored row to contribute.")
  r("trunc()", "categorical", "refused",
    "Refused for the same reason: no order, no CDF, no truncation window.")
  r("trunc()", "cox", "conditional",
    "Runs through the same log-CDF the censoring uses, so trunc(lb = ) is delayed entry. The truncation bound is evaluated against the SAME spline basis the response is, which means a bound outside the boundary knots is clamped to them rather than extrapolated. Untested against an external left-truncated reference.")
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
  r("trials()", "binomial", "works",
    "Also accepted as the glm spelling cbind(successes, failures), which is rewritten to successes | trials(successes + failures); the two fits are identical. The two spellings cannot be combined.")

  ## vint() and vreal() --------------------------------------------------
  # override: "nothing checks this" outranks the permissive blanket
  # defaults (verbose, call_group, double_bar) at the same signature.
  r("vint()", "*", "untested",
    "vint() passes integer covariates to a custom family; only the custom-family path uses them.",
    override = TRUE)
  r("vreal()", "*", "untested",
    "vreal() passes real covariates to a custom family; only the custom-family path uses them.",
    override = TRUE)
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
  # No rescor x kind:method default: every declared method now has an
  # explicit rescor rule, so such a default could never win a pair.
  r("rescor", "fitted", "refused",
    "Refused: fitted() calls single_response() and stops with 'fitted() is not supported yet for multivariate fits'. Predict one response at a time instead: predict(fit, resp = ).")
  r("rescor", "predict", "works",
    "predict() takes resp = to choose the response, and defaults to the first.")
  r("rescor", "simulate", "refused",
    "Refused: simulate() is not supported for multivariate fits yet.")
  r("rescor", "residuals", "refused",
    "Refused: residuals() is not supported for multivariate fits yet.")
  r("rescor", "residuals_osa", "refused",
    "Refused: residuals() is not supported for multivariate fits yet.")
  r("rescor", "emmeans", "refused",
    "Refused: emmeans support is univariate-only for now.")
  # confint() and hypothesis() work on the outer parameter vector,
  # which a multivariate fit has like any other: verified on a
  # two-response gaussian fit, rescor = TRUE included.
  r("rescor", "confint_profile", "works",
    "Verified: the profiled and uniroot intervals for a per-response coefficient agree with the Wald interval on a two-response fit.")
  r("rescor", "hypothesis_profile", "works",
    "Verified: profile likelihood tests address the per-response coefficients by their vcov() names (y1_x and so on).")

  ## mvbf ------------------------------------------------------------------
  r("mvbf", "kind:method", "refused",
    "Refused: the post-fit methods below are univariate-only for now.")
  # residuals_osa x kind:structure claims "untested" at the same
  # signature as the univariate-only refusal above, so the pair is
  # named outright: residuals(type = "osa") goes through the same
  # single_response() guard and stops.
  r("mvbf", "residuals_osa", "refused",
    "Refused: residuals() is not supported for multivariate fits yet, one-step-ahead residuals included.")
  r("mvbf", "fitted", "refused",
    "Refused: fitted() calls single_response() and stops with 'fitted() is not supported yet for multivariate fits'. Predict one response at a time instead: predict(fit, resp = ).")
  r("mvbf", "predict", "works",
    "predict() takes resp = to choose the response, and defaults to the first.")
  # The inference surface reads the outer parameter vector rather than
  # a single response, so the univariate-only refusal above does not
  # reach it. Verified on a two-response gaussian fit.
  r("mvbf", "confint_profile", "works",
    "Verified: the profiled and uniroot intervals for a per-response coefficient agree with the Wald interval on a two-response fit.")
  r("mvbf", "hypothesis_profile", "works",
    "Verified: profile likelihood tests address the per-response coefficients by their vcov() names (y1_x and so on).")
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
  r("|ID|", "gr_cov", "conditional",
    "Works when every term sharing the key is gr(cov = ) over the SAME grouping factor and the SAME matrix: the linked terms merge into one Kronecker block whose covariance is A (x) Sigma, with Sigma unstructured across the merged coefficients. That is the same joint density as the long format, bf(value ~ 0 + trait + (0 + trait | gr(id, cov = A)), sigma ~ 0 + trait), so the two spellings agree to optimizer tolerance. Refused when the key mixes structures (us with gr_cov, cov with prec), spans different grouping specifications, or resolves cov = to different matrices in different formula environments. Put the matrix in data2 so every formula resolves the same object.")
  r("|ID|", "gr_prec", "conditional",
    "Works when every term sharing the key is gr(prec = ) over the SAME grouping factor and the SAME matrix: the merged block precision is Q (x) Sigma^-1, assembled sparsely exactly as for a single correlated-slopes gr(prec = ) term. Refused when the key mixes structures (us with gr_prec, cov with prec), spans different grouping specifications, or resolves prec = to different matrices in different formula environments. Put the matrix in data2 so every formula resolves the same object.")

  ## Student-t latents, gr(dist = "student") --------------------------------
  r("|ID|", "group:student_blocks", "refused",
    "Refused: a merged |ID| block is assembled as a gaussian one. Write the merged coefficients as a single term instead: (x1 + x2 | gr(g, dist = \"student\")) is the same multivariate-t block, since the t's mixing variable is per level and already shared across the level's coefficients.")
  r("mm()", "group:student_blocks", "refused",
    "Refused: an mm() row loads several levels at once, and the t's mixing variable is one per level, so the row has no single value of it. brms accepts mm(..., dist = \"student\") as a hierarchical construction it samples; there is no closed-form marginal density here to hand the Laplace machinery.")
  r("group:student_blocks", "gr_cov", "refused",
    "Refused: gr(g, cov = A, dist = \"student\") correlates the LEVELS through A while the t's mixing variable is per level, so the joint density over the field is not a multivariate t. brms writes dfm .* (sd * (Lcov * z)) for this and samples it; it has no closed form to marginalize.")
  r("group:student_blocks", "gr_prec", "refused",
    "Refused for the same reason as gr(cov = ): a precision matrix over the levels and a per-level mixing variable do not compose into a density the Laplace approximation can be taken over.")
  r("group:student_blocks", "group:post_fit", "works",
    "Verified: the Laplace machinery is distribution-agnostic, so ranef(), its conditional variances and sdreport() standard errors need no change. simulate() and frm_simulate() draw a multivariate t with one chi-square per level. predict(allow_new_levels = TRUE) inflates the unseen level's variance by nu/(nu-2); the interval is still built as a gaussian one, so it carries the right variance but not the right far-tail quantile.")
  r("quadrature", "group:student_blocks", "conditional",
    "Allowed for one-dimensional blocks, and RECOMMENDED there: the Gauss-Kronrod rule marginalizes a scalar t latent EXACTLY, where the Laplace default is approximate. Verified against adaptive Gauss-Hermite quadrature to 1e-6 in the log-likelihood and in every estimate. Correlated slopes are refused, as they are for a gaussian block.")
  r("REML", "group:student_blocks", "works",
    "Verified against a two-stage exact quadrature: REML stacks one more Laplace dimension on the t latent, and the objective error roughly doubles while staying small (0.11 against 0.05 at nu = 3 over 15 groups). The variance component moves in the direction REML exists to move it.")

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
  r("mixture", "kind:covstruct", "untested", "")
  r("mixture_mvn", "simulate", "works",
    "Since v0.36 a structured simulator (fam$sim_ctx) draws a class per row from the gating weights and then a D-variate normal about that class's mean with that class's covariance, assembled from the family-level extras by the same sigma() the likelihood uses. A draw is an n x D matrix. simulate(), posterior_predict() and frm_simulate() all reach the same implementation.")
  r("mixture_mvn", "mvbf", "refused",
    "Refused: the family already takes a matrix response.")
  r("mixture_mvn", "kind:aterm", "untested", "")

  ## specials ------------------------------------------------------------------
  r("mo()", "kind:family", "conditional",
    "The monotonic variable must be an ordered factor or an integer with at least 3 categories. An interaction multiplier (mo(x):z) must be a single numeric column: factor and character multipliers are refused, because the simplex carries one coefficient and a contrast expansion has no column to go in. Expand the factor to numeric indicators first. brms refuses the same spelling (brms#1828).")
  # A rule may not name the same feature twice: the resolved table
  # holds unordered pairs of DISTINCT features, so a self-pair row
  # could never match. mo() x mo() lives in the note below instead.
  r("mo()", "mi_pred()", "refused",
    "Refused: mo() cannot interact with another mo() or mi() term. Two mo() terms may appear in one formula, but not crossed with each other.")
  # override: the sparse_x / autoscale sweeps claim "works" for every
  # special at the same signature. Nothing checks mo(), s() or gp()
  # under those modes, so the untested claim governs.
  r("mo()", "kind:mode", "untested", "", override = TRUE)
  r("mo()", "REML", "works", "Verified by a tiny fit.")
  r("mo()", "profile", "works", "Verified by a tiny fit.")
  r("mo()", "predict", "conditional",
    "New data must stay inside the fitted category range; unknown categories are refused.")
  r("cs_pred()", "kind:family", "refused",
    "Refused: cs() needs an sratio, cratio, or acat family.")
  r("cs_pred()", "group:ordinal_cs", "works", "")
  r("cs_pred()", "cumulative", "refused",
    "Refused: category-specific effects are not identified under the cumulative parameterization.")
  r("gp_pred()", "kind:mode", "untested", "", override = TRUE)
  # gp()'s own arity limits used to sit on a gp_pred() x gp_pred()
  # self-pair row, which the resolved table cannot hold; every gp()
  # model has a family, so they are stated here instead.
  r("gp_pred()", "kind:family", "works",
    "gp() takes 1 to 3 variables. The arguments k, c, and iso are evaluated in the formula environment, and c may be a vector with one entry per dimension. Several gp() terms may appear in one formula.")
  r("s()", "kind:mode", "untested", "", override = TRUE)
  r("s()", "kind:family", "works", "")
  r("t2()", "kind:family", "works", "")

  ## covariance structures -------------------------------------------------------
  #
  # Each structure states its own condition against everything. The
  # condition is a property of the structure, not of what it is paired
  # with, so these carry override = TRUE: at the blanket signature they
  # must outrank the permissive and untested defaults above, and a
  # structure that needs num_factor() coordinates or a level order must
  # not be reported as unconditionally working. They are written one
  # name per rule rather than by group, so that every structure's
  # condition carries the same authority.
  r("group:positional", "kind:special", "conditional",
    "Positional structures take the lag from the level order, not from the level value.")
  r("ar1", "*", "conditional",
    "ar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns. Use ou() over num_factor() for irregular spacing.",
    override = TRUE)
  r("hetar1", "*", "conditional",
    "hetar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns.",
    override = TRUE)
  r("ou", "*", "conditional",
    "ou() is the irregular-spacing structure. Build the levels with num_factor() so the metric distance is recoverable.",
    override = TRUE)
  for (cs in frmtmb_compat_groups_lst$spatial) {
    r(cs, "*", "conditional",
      "Spatial structures need coordinates built with num_factor(x, y).",
      override = TRUE)
  }
  r("gr_prec", "*", "conditional",
    "gr(prec = Q) takes correlated slopes; the block precision is the Kronecker product of Q and the inverse term covariance, so it stays as sparse as Q. Q needs dimnames covering every grouping level, and belongs in data2 = list(Q = Q). Terms sharing an |ID| key over the same factor and the same Q merge into one such block.",
    override = TRUE)
  r("car", "*", "conditional",
    "car(M, gr = g, type = ) is a predictor special, not a bar term. M is a symmetric adjacency matrix with dimnames (rownames, colnames, or both, which then have to agree) covering every location; entries must be present and non-negative, and non-zero weights are binarized. type = \"escar\" is the proper CAR, \"icar\"/\"esicar\" the intrinsic one under a soft sum-to-zero constraint (con_sd), \"bym2\" the scaled mixture; escar needs every location to have a neighbor. M belongs in data2 = list(M = M).",
    override = TRUE)
  r("spde", "*", "conditional",
    "spde(fem, gr = node) is a predictor special taking a mesh's finite-element matrices (M0/M1/M2 or c0/g1/g2) as fixed data; gr maps observations onto mesh nodes BY ROW NUMBER (whole numbers in 1..nrow(M0), as integers or as a factor/character spelling of them), because the matrices carry no dimnames to match labels against. Unobserved nodes keep their column; a general projector matrix is not supported yet. The matrices belong in data2 = list(fem = fem).",
    override = TRUE)
  for (cs in frmtmb_compat_groups_lst$gmrf) {
    r(cs, "simulate", "works",
      "Draws come from the block's own fitted precision.")
    r(cs, "predict", "conditional",
      "Locations outside the fitted set have no marginal variance of their own; allow_new_levels predicts them at the population level with no block contribution.")
  }
  r("gr_cov", "*", "conditional",
    "gr(cov = A) accepts correlated slopes; the block covariance is the Kronecker product of A and the term covariance. A needs dimnames covering every grouping level, and belongs in data2 = list(A = A). Terms sharing an |ID| key over the same factor and the same A merge into one such block, which is the same model as writing the traits long with a single gr() term.",
    override = TRUE)
  r("equalto", "*", "conditional",
    "equalto(x + 0 | g, V) fixes the term covariance to V, which must be square and match the term dimension, and belongs in data2 = list(V = V).",
    override = TRUE)
  r("rr", "*", "conditional",
    "rr() gives a reduced-rank block; the rank d must not exceed the term dimension.",
    override = TRUE)
  r("rr", "sparse_x", "works",
    "Verified: identical estimates with and without a sparse fixed-effect design.")
  r("rr", "REML", "works", "Verified by a tiny fit.")
  r("smooth", "*", "conditional",
    "smooth is the internal structure behind s() and t2(); it is not written directly in a formula.",
    override = TRUE)
  for (cs in frmtmb_compat_groups_lst$latent_gp) {
    r(cs, "*", "conditional",
      "gp() and hsgp() are predictor specials, not bar terms. Write gp(x), not (gp(x) | g).",
      override = TRUE)
  }

  ## post-fit methods --------------------------------------------------------------
  r("simulate", "kind:family", "works",
    "The family supplies a simulator.")
  r("simulate", "group:no_simulator", "refused",
    "Refused: this family has no simulator yet.")
  r("simulate", "group:ordinal", "works",
    "Draws come back as an ordered factor carrying the response's own levels, not as 1..K codes.")
  r("simulate", "multinomial", "works",
    "Draws come back as an n x K count matrix carrying the response's column names; needs trials().")
  r("simulate", "categorical", "works",
    "Draws come back as an UNORDERED factor carrying the response's own levels: the categories have no order, so restoring them as an ordered factor would claim one the model never used.")
  r("simulate", "von_mises", "works",
    "Best and Fisher's rejection sampler, vectorized over per-row mu and kappa, so a distributional kappa simulates. RTMBdist::rvm() takes scalar parameters only and is not used.")
  r("simulate", "cox", "refused",
    "Refused: drawing a survival time means inverting the cumulative baseline hazard, which this family does not carry a quantile function for. simulate(), posterior_predict() and frm_simulate() each say so in their own words and then repeat the family's reason.")
  r("residuals_osa", "kind:family", "conditional",
    "One-step-ahead residuals need the family to register its observation through OBS().")
  r("residuals_osa", "categorical", "refused",
    "Refused with residuals() as a whole: a one-step-ahead residual is a CDF value, and a nominal response has no CDF.")
  r("residuals", "categorical", "refused",
    "Refused: the categories carry no order, so no residual has a scale to live on. Compare fitted(fit), the n x K category probabilities, against the observed categories instead.")
  r("residuals_osa", "von_mises", "refused",
    "Refused upstream: RTMBdist::dvm() rejects the osa observation object, because a wrapped support has no one-step CDF on the line.")
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
  r("emmeans", "group:ordinal", "conditional",
    "Works on the LATENT linear predictor, emmeans's mode = \"latent\" convention for clm-like models: the intercept is dropped there (the K-1 thresholds take its place), so contrasts are on the latent scale and absolute means carry no threshold offset. For category probabilities use predict(fit, type = \"response\") or conditional_effects(), which are on a different scale from these means.")
  r("confint_profile", "kind:mode", "untested", "")
  r("hypothesis_profile", "kind:mode", "untested", "")
  r("predict", "kind:family", "conditional",
    "Rank-deficient designs drop aliased columns at fit time. New data that is not estimable from the retained columns predicts NA and warns.")
  r("predict", "group:ordinal", "conditional",
    "type = \"response\" returns an n x K matrix of category probabilities (rows summing to 1, columns named by the response's own levels), not a vector: an ordinal response has no mean. It equals fitted(). cs() terms are honored and re-evaluated on newdata. type = \"link\" gives the latent predictor, which is where se.fit is available; se.fit is refused on the response scale.")
  r("predict", "categorical", "conditional",
    "type = \"response\" returns an n x K matrix of category probabilities, columns named by the response's own levels and rows summing to 1, exactly as for the ordinal families; it equals fitted(). se.fit is refused there. Each category's latent predictor is predict(type = \"link\", dpar = \"mu<Level>\"), which is where se.fit works.")
  r("predict", "cox", "conditional",
    "type = \"response\" and fitted() are refused: a survival time has no mean the censored rows identify. type = \"link\" gives the log hazard ratio, and cox_baseline() the fitted baseline weights.")
  r("fitted", "kind:family", "conditional",
    "Needs a family with a mean function.")
  r("fitted", "categorical", "conditional",
    "Returns the n x K matrix of category probabilities, not a vector: a nominal response has no mean, so the modelled response is the category distribution. The predict(type = \"response\") == fitted() identity holds.")
  r("fitted", "cox", "refused",
    "Refused: a survival time has no mean on the response scale here. Use predict(type = \"link\") for the log hazard ratio.")
  r("fitted", "von_mises", "conditional",
    "Returns the mean DIRECTION in radians on (-pi, pi], which is what brms's posterior_epred() reports for this family; a circular response has no arithmetic mean.")
  r("fitted", "group:ordinal", "conditional",
    "Returns the same n x K matrix of category probabilities predict(type = \"response\") returns, not a vector: an ordinal response has no mean, so the modelled response is the category distribution. The predict(type = \"response\") == fitted() identity holds. The latent linear predictor is predict(fit, type = \"link\"), which is also what emmeans and insight see.")
  r("residuals", "group:ordinal", "conditional",
    "\"response\" and \"pearson\" score the categories by the same codes 1..K the likelihood uses: y - sum_k k * P(y = k), standardized by that distribution's own sd. That is a residual on a SCORE, not on the ordinal scale; \"osa\" and dharma_residuals() use only the order. \"deviance\" is refused, as for every family without a standard unit deviance.")

  ## R-side residual correlation -----------------------------------------------
  # Everything refused here is refused for one reason: the likelihood
  # is a joint density over each group, so it no longer factorizes into
  # per-row contributions. brms refuses the same core set.
  r("group:autocor", "kind:family", "refused",
    "Refused: a residual correlation needs a family with a real residual. brms accepts the same spelling for other families but fits a different model there - a latent gaussian AR process added to the linear predictor - which is spelled here as a random effect over the time factor: + ar1(factor(week) + 0 | subj), or toep()/us() for a freer lag structure.")
  r("group:autocor", "group:autocor_families", "works",
    "Written as a formula term - ar(week, subj, cov = TRUE), cosy(gr = subj), unstr(week, subj) - it makes the residuals of one group a single correlated draw, y_g ~ N(mu_g, D R D) with D the diagonal of that group's sigma values; student() gets the multivariate-t analog. These are exactly the two families brms treats this way. brms's default cov = FALSE (the residual-regression formulation) is a different likelihood and is refused; the call must say cov = TRUE. The lag is the distance between the rows' positions in the GLOBAL set of time levels, so a group missing a time point gets the wider lag (nlme's reading, not brms's). Validated against nlme::gls (corAR1, corARMA, corCompSymm, corSymm) under ML and REML: log-likelihoods agree to 1e-9 or better and the correlation parameters to 1e-5 or better.")
  r("group:autocor", "student", "conditional",
    "The multivariate-t has one shape parameter per group, so nu must be constant; a predicted nu ~ ... is refused. The density is brms's multi_student_t with scale matrix D R D, verified against mvtnorm::dmvt exactly.")
  r("group:autocor", "kind:aterm", "refused",
    "Refused: the group's density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with 'Invalid addition arguments for this model'.")
  r("group:autocor", "rescor", "refused",
    "Refused: both describe the residual covariance - one across time, one across responses - and the joint structure is their Kronecker product, which is not implemented. brms refuses the same pair.")
  r("group:autocor", "mixture", "refused",
    "Refused: a mixture likelihood has no single residual to correlate. The term is rejected as sitting on mu1 rather than mu, which is also how brms rejects it.")
  r("group:autocor", "mixture_mvn", "refused",
    "Refused for the same reason as mixture().")
  r("group:autocor", "nl", "refused",
    "Refused: a nonlinear mu is arbitrary R code, so the term would be evaluated rather than read. brms reaches this model through acformula(), which has no analog here.")
  r("group:autocor", "quadrature", "refused",
    "Refused: the Gauss-Kronrod rule integrates a random effect against per-observation densities, and this residual is a joint density over each group.")
  r("group:autocor", "REML", "works",
    "The correlation parameters are covariance parameters and stay in the outer problem, exactly as theta does. Verified against gls(method = \"REML\") and lme(method = \"REML\") to 1e-10 in the log-likelihood.")
  r("group:autocor", "profile", "works",
    "Verified: the same log-likelihood as the plain ML fit.")
  r("group:autocor", "sparse_x", "works", "Verified by a tiny fit.")
  r("group:autocor", "autoscale", "works", "Verified by a tiny fit.")
  r("group:autocor", "bounds", "works",
    "The parameters are the thetaac_* rows of the outer vector and can be bounded by that name.")
  r("group:autocor", "prior", "conditional",
    "Priors on the fixed effects and on random-effect covariance parameters work as usual. set_prior() cannot target the residual-correlation parameters themselves yet; bounds on thetaac_* are the available lever.")
  r("group:autocor", "kind:covstruct", "works",
    "Random effects alongside a correlated residual are the point of the feature: the marginal likelihood is a Laplace approximation over the modes with the multivariate residual density inside. Verified against nlme::lme(random = ~ 1 | subj, correlation = corAR1()) under ML and REML.")
  r("group:autocor", "mvbf", "works",
    "One residual correlation term per response, each with its own parameters. Only with rescor = FALSE.")
  r("group:autocor", "|ID|", "untested", "")
  r("group:autocor", "simulate", "works",
    "simulate() draws one correlated residual per group (a Cholesky factor of that group's correlation submatrix applied to standard normal, or scaled-t, innovations), so the draws carry the fitted autocorrelation. dharma_residuals() therefore works too. Since v0.36 the residual draw is part of the simulator contract rather than a branch inside simulate(), so posterior_predict() and frm_simulate() draw it as well; frm_simulate() needs thetaac in the internal newparams spelling.")
  r("group:autocor", "residuals_osa", "refused",
    "Refused: one-step-ahead residuals need the taped density of one observation given the previous ones, and the tape holds a joint density per group. Use type = \"pearson\", which divides by the marginal residual SD (the diagonal of the residual covariance is sigma^2, because R is unit-diagonal), or dharma_residuals().")
  r("group:autocor", "kind:method", "untested", "")
  r("group:autocor", "fitted", "works",
    "The mean structure is untouched: the term changes the residual density, not the linear predictor.")
  r("group:autocor", "predict", "works",
    "Unchanged, newdata and se.fit included: se.fit is the uncertainty of the MEAN, which the residual correlation does not enter.")
  r("group:autocor", "residuals", "conditional",
    "\"response\" and \"pearson\" are unchanged, and pearson is still the right standardization: R is unit-diagonal, so the marginal residual SD is sigma. They do NOT decorrelate the residuals, so a plot of them against time still shows the fitted autocorrelation; that is the intended reading.")
  r("group:autocor", "confint_profile", "works",
    "tmbprofile() addresses the residual-correlation parameters under their thetaac_* names.")
  r("group:autocor", "kind:special", "untested", "")
  r("group:autocor", "kind:autocor", "refused",
    "Refused: a response has one residual covariance, so it carries one such term. brms refuses the same with 'Can only model one time-series term'.")

  ## formula grammar --------------------------------------------------------------
  # The two permissive grammar defaults are declared with the other
  # weak blanket rules near the top of the table; only the refusals
  # belong here, at the authoritative end of the order.
  r("bar_crossing", "*", "refused",
    "Refused: a bar term crossed with * or : (as in x * (1 | g)) is not a random-effect specification (lme4#196). Write the crossing inside the bar: (x | g). This spelling was once accepted with the crossing silently dropped.",
    override = TRUE)
  r("double_bar", "diag", "works",
    "diag is the structure (f || g) resolves to.")
  r("bar_crossing", "kind:covstruct", "refused",
    "Refused for every covariance structure. The crossing must go inside the bar.")
  # The other two grammar rows are "*" rules of the same shape as the
  # bar_crossing one, so without this the crossing refusal and the
  # call_group / double_bar support rules would tie and file order
  # would decide. The parser refuses the crossing before it ever looks
  # at how the group or the bar is spelled: verified for
  # x * (1 | factor(g)) and x * (1 + x || g), both of which stop with
  # "A random-effect term cannot be crossed with '*'".
  r("bar_crossing", "kind:grammar", "refused",
    "Refused whatever the bar contains: the crossing is rejected before the grouping factor or the double bar is read. Write the crossing inside the bar.")

  ## multi-membership ------------------------------------------------------
  # mm() changes the Z MATRIX and nothing else: the block is an
  # ordinary us/diag block over the pooled level set, so the family,
  # the addition terms and every post-fit method meet the same object
  # they meet for (1 | g). What mm() cannot do is share a block with
  # anything that indexes ONE level per observation row.
  r("mm()", "kind:family", "works",
    "The membership design is built before the family sees it, exactly as an ordinary grouping factor's is.")
  r("mm()", "kind:aterm", "works",
    "Addition terms change the likelihood; the membership design changes the predictor.")
  r("mm()", "kind:special", "works",
    "Smooths, Gaussian processes and monotonic terms are separate additive terms and separate blocks.")
  r("mm()", "group:post_fit", "works",
    "Verified: an mm block is an ordinary block, so ranef(), VarCorr(), ngrps(), fitted(), simulate() and residuals() read it with no multi-membership branch. predict() on newdata rebuilds the weighted rows; a membership level that is new needs allow_new_levels = TRUE, and the row's remaining members still contribute their fitted effects.")
  r("mm()", "kind:covstruct", "refused",
    "Refused: mm() supports the default (us) and diag structures only. Every other structure describes a covariance over the block's levels, and the pooled membership levels have no ordering, no coordinates and no relationship matrix for one to be defined on. gr(mm(...), cov = ) is refused for the same reason.")
  r("mm()", "us", "works",
    "The default: one unstructured covariance over the term's coefficients, shared by every pooled membership level.")
  r("mm()", "diag", "works",
    "diag(x | mm(g1, g2)) is the spelling of brms's mm(..., cor = FALSE); (x || mm(g1, g2)) expands to the same thing.")
  r("mm()", "|ID|", "refused",
    "Refused: a merged |ID| block indexes one level set per observation row, and an mm() row loads several levels at once.")
  r("mm()", "REML", "works",
    "Verified: REML integrates the fixed effects and leaves the block alone.")
  r("mm()", "quadrature", "untested",
    "An (1 | mm(g1, g2)) block passes the scalar-intercept guard, but nothing checks that the Gauss-Kronrod rule marginalizes a design whose rows load several levels at once. Use the Laplace default.")
  r("mm()", "mvbf", "works",
    "Each response builds its own membership design; the blocks do not interact.")
  r("mm()", "mmc()", "works",
    "mmc(x1, x2) is ONE random-slope coefficient of the mm block whose covariate value is member specific: member k takes argument k. Verified against brms's Z_..._k arrays.")
  r("mmc()", "*", "conditional",
    "mmc() only means something on the left of a multi-membership bar, where it supplies one covariate value per member. Anywhere else it is refused, including over a single-membership grouping factor. Inside an mm() term it composes like any other random-slope column.",
    override = TRUE)

  # A covariance structure's own conditions (num_factor() coordinates
  # for exp/gau/mat, level order for ar1, a rank for rr) hold whatever
  # the response distribution is, so there is deliberately NO
  # <covstruct> x kind:family rule granting "works" here. Such a rule
  # would outrank each structure's own "<name> x *" condition and
  # replace it with an unconditional claim.

  out <- do.call(rbind, c(rows,
                          lapply(frmtmb_compat_contrib$rules,
                                 function(mk) mk())))
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------- resolution

frmtmb_compat_statuses <- c("works", "conditional", "refused",
                            "broken", "untested")

#' Specificity of one pattern. Higher beats lower.
#'
#' @noRd
compat_spec <- function(pat) {
  ifelse(pat == "*", 0L,
         ifelse(startsWith(pat, "kind:"), 1L,
                ifelse(startsWith(pat, "group:"), 2L, 3L)))
}

#' A rule's precedence, as one comparable integer: the two sides'
#' specificities sorted descending, then read as a two-digit base-4
#' number. Lexicographic order on the sorted pair is exactly numeric
#' order on the result, so the resolver can compare with `>=`. See the
#' PRECEDENCE note above `frmtmb_compat_rules_tbl()`.
#'
#' @noRd
compat_sig_rank <- function(pa, pb) {
  a <- compat_spec(pa)
  b <- compat_spec(pb)
  4L * pmax(a, b) + pmin(a, b)
}

#' Does `pat` cover the feature at position i of the feature table?
#'
#' @noRd
compat_match <- function(pat, name, kind) {
  if (pat == "*") return(rep(TRUE, length(name)))
  if (startsWith(pat, "kind:")) return(kind == substring(pat, 6L))
  if (startsWith(pat, "group:")) {
    g <- frmtmb_compat_groups_lst[[substring(pat, 7L)]]
    return(name %in% g)
  }
  name == pat
}

#' A rule is unordered: it may match a pair either way round.
#'
#' @noRd
compat_hit <- function(pa, pb, pairs) {
  (compat_match(pa, pairs$feature_a, pairs$kind_a) &
     compat_match(pb, pairs$feature_b, pairs$kind_b)) |
    (compat_match(pa, pairs$feature_b, pairs$kind_b) &
       compat_match(pb, pairs$feature_a, pairs$kind_a))
}

#' Unordered pairs of distinct features, minus family x family: a model
#' carries exactly one family, so those pairs mean nothing.
#'
#' @noRd
frmtmb_compat_pairs_tbl <- function() {
  ft <- frmtmb_compat_features_tbl()
  ij <- utils::combn(nrow(ft), 2L)
  keep <- !(ft$kind[ij[1, ]] == "family" & ft$kind[ij[2, ]] == "family")
  ij <- ij[, keep, drop = FALSE]
  data.frame(
    feature_a = ft$name[ij[1, ]], kind_a = ft$kind[ij[1, ]],
    feature_b = ft$name[ij[2, ]], kind_b = ft$kind[ij[2, ]],
    stringsAsFactors = FALSE)
}

#' Which rule wins each pair. Resolution always runs over the whole pair
#' table, never over the slice a caller asked for, so a rule's
#' precedence is judged against the registry rather than against the
#' query.
#'
#' @noRd
compat_resolve <- function(pairs, rules) {
  rank <- compat_sig_rank(rules$feature_a, rules$feature_b)
  best <- rep(-1L, nrow(pairs))
  win <- rep(NA_integer_, nrow(pairs))
  for (k in seq_len(nrow(rules))) {
    # >= so a later rule of equal precedence overrides an earlier one
    hit <- compat_hit(rules$feature_a[k], rules$feature_b[k], pairs) &
      rank[k] >= best
    win[hit] <- k
    best[hit] <- rank[k]
  }
  list(win = win, best = best, rank = rank)
}

#' Pairs whose winning precedence is shared by rules that disagree about
#' the status. Such a pair is decided by file order, which is the defect
#' the sorted-pair precedence exists to remove; the registry must have
#' none of them unless the winning rule says it is a deliberate
#' override. Run from test-compat.R rather than at load time: it costs a
#' pass over every rule x pair combination.
#'
#' @noRd
frmtmb_compat_validate <- function() {
  pairs <- frmtmb_compat_pairs_tbl()
  rules <- frmtmb_compat_rules_tbl()
  res <- compat_resolve(pairs, rules)
  ovr <- if (is.null(rules$override)) rep(FALSE, nrow(rules)) else
    rules$override
  seen <- integer(nrow(pairs))       # rule index of the last top hitter
  bad <- integer(0)
  for (k in seq_len(nrow(rules))) {
    hit <- compat_hit(rules$feature_a[k], rules$feature_b[k], pairs) &
      res$rank[k] == res$best
    if (!any(hit)) next
    clash <- which(hit & seen > 0L)
    if (length(clash)) {
      clash <- clash[rules$status[seen[clash]] != rules$status[k]]
      bad <- c(bad, clash)
    }
    seen[hit] <- k
  }
  bad <- unique(bad)
  bad <- bad[!ovr[res$win[bad]]]
  w <- res$win[bad]
  data.frame(
    feature_a = pairs$feature_a[bad], feature_b = pairs$feature_b[bad],
    status = rules$status[w],
    # paste() would recycle a zero-length index up to length one
    rule = if (length(w)) paste(rules$feature_a[w], "x",
                                rules$feature_b[w]) else character(0),
    stringsAsFactors = FALSE)
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
#' feature (specificity 0), `"kind:<kind>"` for a whole kind (1),
#' `"group:<group>"` for a named feature set (2), or a bare feature
#' name (3).
#'
#' A pair takes the status of the rule with the highest precedence,
#' where precedence is the two sides' specificities sorted descending
#' and compared lexicographically: `(3,3)` beats `(3,2)` beats `(3,1)`
#' beats `(3,0)` beats `(2,2)` beats `(2,1)`, and so on. A rule that is
#' strictly more specific on one side and no less specific on the other
#' therefore always wins. Rules can tie only when their two
#' specificities match exactly; the later rule wins those, so an
#' override is appended rather than inserted, and a tie between rules
#' that disagree about the status is a registry defect the test suite
#' rejects.
#'
#' Use [frm_compat()] to read the resolved answer for a pair. Use this
#' function to see which rule is doing the work.
#'
#' @return A data frame with columns `feature_a`, `feature_b`,
#'   `status`, `note`, and `override`.
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
#'   pair involving that feature, or neither for the whole table. Both
#'   accept a vector, which gives every pair in the cross of the two
#'   sides; an empty vector is an error rather than an empty answer.
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
  # An empty selector is a caller mistake, not a query for nothing: it
  # used to fall through the recycling below and return an empty table,
  # which reads exactly like "this feature is involved in no pairs".
  check_features <- function(x, arg) {
    if (is.null(x)) return(NULL)
    if (!is.character(x)) {
      stop(arg, " must be a character vector of feature names.",
           call. = FALSE)
    }
    if (!length(x)) {
      stop(arg, " is empty. Supply at least one feature name, or NULL ",
           "for every feature.", call. = FALSE)
    }
    bad <- setdiff(x, ft$name)
    if (length(bad)) {
      stop("Unknown feature: '", bad[1], "'. See frm_compat_features().",
           call. = FALSE)
    }
    x
  }
  feature_a <- check_features(feature_a, "feature_a")
  feature_b <- check_features(feature_b, "feature_b")
  if (!is.null(status)) {
    bad <- setdiff(status, frmtmb_compat_statuses)
    if (length(bad)) {
      stop("Unknown status: ", paste(bad, collapse = ", "),
           ". Statuses are: ",
           paste(frmtmb_compat_statuses, collapse = ", "),
           call. = FALSE)
    }
  }

  # Resolve before subsetting so a rule's precedence is judged against
  # the whole table, never against the slice the caller asked for.
  pairs <- frmtmb_compat_pairs_tbl()
  rules <- frmtmb_compat_rules_tbl()
  win <- compat_resolve(pairs, rules)$win
  pairs$status <- rules$status[win]
  pairs$note <- rules$note[win]

  # Both sides take a vector, and the answer is the full cross: every
  # pair with one member in feature_a and the other in feature_b.
  sel <- rep(TRUE, nrow(pairs))
  if (!is.null(feature_a) && !is.null(feature_b)) {
    sel <- (pairs$feature_a %in% feature_a &
              pairs$feature_b %in% feature_b) |
      (pairs$feature_a %in% feature_b & pairs$feature_b %in% feature_a)
  } else {
    one <- feature_a %||% feature_b
    if (!is.null(one)) {
      sel <- pairs$feature_a %in% one | pairs$feature_b %in% one
    }
  }
  if (!is.null(status)) sel <- sel & pairs$status %in% status
  out <- pairs[sel, , drop = FALSE]

  # read the queried side first, so a one-feature slice reads down a
  # single column
  lead <- feature_a %||% feature_b
  if (!is.null(lead)) {
    flip <- out$feature_b %in% lead & !(out$feature_a %in% lead)
    if (any(flip)) {
      out[flip, c("feature_a", "kind_a", "feature_b", "kind_b")] <-
        out[flip, c("feature_b", "kind_b", "feature_a", "kind_a")]
    }
  }
  rownames(out) <- NULL
  out
}
