# Where each brms spelling lives (dev/core-boundary.md, guardrail 1).
#
# brms parity is measured on core plus the extensions, so a brms name
# may resolve in core, resolve in a named extension, or be refused by
# name with the reason. What must never happen is the fourth outcome: a
# brms spelling that reaches a generic "not supported" list, or "could
# not find function", which reads as "frmtmb cannot do this" when one
# package away it can. Families have extension_family_homes
# (R/families.R); this file holds the same tables for formula specials
# and addition terms, and the brms vocabulary that test-guardrail1.R
# checks each name of.

#' brms 2.23.0's formula vocabulary: the families `brmsfamily()`
#' accepts (its `.family_*` functions, less the `.family_info()`
#' helper), the addition terms (its exported `resp_*()` functions), and
#' the predictor and group-level specials it exports. Fixed data, so
#' that the guardrail test runs without brms; `test-guardrail1.R`
#' compares it with an installed brms.
#'
#' @noRd
brms_vocabulary <- list(
  families = c(
    "acat", "asym_laplace", "bernoulli", "beta", "beta_binomial",
    "binomial", "categorical", "com_poisson", "cox", "cratio",
    "cumulative", "custom", "dirichlet", "dirichlet_multinomial",
    "dirichlet2", "discrete_weibull", "exgaussian", "exponential",
    "frechet", "gamma", "gaussian", "gen_extreme_value", "geometric",
    "hurdle_cumulative", "hurdle_gamma", "hurdle_lognormal",
    "hurdle_negbinomial", "hurdle_poisson", "inverse.gaussian",
    "logistic_normal", "lognormal", "multinomial", "negbinomial",
    "negbinomial2", "poisson", "shifted_lognormal", "skew_normal",
    "sratio", "student", "von_mises", "weibull", "wiener", "xbeta",
    "zero_inflated_asym_laplace", "zero_inflated_beta",
    "zero_inflated_beta_binomial", "zero_inflated_binomial",
    "zero_inflated_negbinomial", "zero_inflated_poisson",
    "zero_one_inflated_beta"),
  aterms = c("bhaz", "cat", "cens", "dec", "index", "mi", "rate", "se",
             "subset", "thres", "trials", "trunc", "vint", "vreal",
             "weights"),
  specials = c("s", "t2", "gp", "mo", "me", "mi", "cs", "cse", "car",
               "sar", "fcor", "ar", "ma", "arma", "cosy", "unstr", "gr",
               "mm", "mmc")
)

# Formula specials an extension supplies rather than core, keyed by the
# name a formula writes. One entry per line, so that entries added by
# different changes merge.
extension_special_homes <- c(
  car = "frmtmb.spatial",
  spde = "frmtmb.spatial"
)

# Addition terms an extension supplies rather than core.
extension_aterm_homes <- c(
  dec = "frmtmb.eam"
)

# brms families frmtmb refuses by name, with the frmtmb spelling;
# family_ctor() raises it.
brms_family_refusals <- c(
  custom = paste0(
    "brms's custom family is built with custom_family(), not with ",
    "brmsfamily(\"custom\"); frmtmb has the same constructor: pass ",
    "family = custom_family(...)")
)

# brms formula specials frmtmb refuses by name. parse_linpred() raises
# the refusal with the term and the residual correlation terms there are.
brms_special_refusals <- c(
  sar = "is a brms autocorrelation term that frmtmb does not support",
  fcor = "is a brms autocorrelation term that frmtmb does not support"
)

# Where the bhaz() refusal sends a user: the family that takes brms's
# baseline-hazard basis. The one place to change when cox() moves to an
# extension package (the refusal reads it).
bhaz_family_pointer <- paste0(
  "family = cox(df = 5, degree = 3, intercept = TRUE), which is brms's ",
  "bhaz() default")

# brms addition terms frmtmb refuses, with the reason and the frmtmb
# spelling where there is one.
brms_aterm_refusals <- c(
  bhaz = paste0(
    "it sets brms's M-spline basis for the baseline hazard of a cox() ",
    "model, and frmtmb sets the same basis in the family: ",
    bhaz_family_pointer, ". A group-specific baseline, bhaz(gr = ), has ",
    "no frmtmb equivalent")
)

#' The extension that supplies a formula special, or `NULL`.
#'
#' @noRd
extension_special_home <- function(name) {
  if (!name %in% names(extension_special_homes)) return(NULL)
  unname(extension_special_homes[[name]])
}

#' The refusal for a spelling whose home extension has not registered
#' it, in the words the family refusal uses (`family_ctor()`).
#'
#' @noRd
stop_extension_home <- function(name, home, what) {
  frm_stop("The ", name, "() ", what, " is in the extension package ",
           home, ", which is not loaded. Install it and run library(",
           home, "), then use ", name, "()", call. = FALSE)
}

#' Refuse a formula term that calls a special an extension supplies,
#' when that extension has not registered it. A registered special is
#' a covariance structure's field term, so it is in the registry when
#' its package is loaded; the call position is read with
#' `called_functions()`, so a data column of the same name is not
#' mistaken for the special.
#'
#' @noRd
check_extension_specials <- function(tm) {
  nms <- setdiff(names(extension_special_homes), cs_names())
  if (!length(nms)) return(invisible(NULL))
  hit <- called_functions(tm, nms)
  if (length(hit)) {
    stop_extension_home(hit[[1L]], extension_special_home(hit[[1L]]),
                        "special")
  }
  invisible(NULL)
}

#' Refuse an addition term that is not core's and not registered, by
#' name: the extension that supplies it, or why frmtmb refuses it.
#' `NULL` when neither table names it, so the caller's general refusal
#' runs.
#'
#' @noRd
check_unsupported_aterm <- function(nm) {
  if (nm %in% names(extension_aterm_homes)) {
    stop_extension_home(nm, unname(extension_aterm_homes[[nm]]),
                        "addition term")
  }
  if (nm %in% names(brms_aterm_refusals)) {
    frm_stop("Addition term `", nm, "()` is not supported: ",
             brms_aterm_refusals[[nm]], call. = FALSE)
  }
  invisible(NULL)
}
