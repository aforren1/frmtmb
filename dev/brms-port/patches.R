# Documented spelling changes, applied only in the second pass.
#
# Two kinds. AUTO_RETRY rules are generic: they look at the error frm()
# raised and rewrite the call the way the porting guide says to. PATCH
# entries are per-expression edits that no rule could infer (start values
# for a nonlinear fit, a formula the grammar spells differently).
#
# Nothing here fixes package code. Each entry is a claim that a user with
# the migration guide in hand could have made the same edit.

add_arg <- function(src, arg, value) {
  e <- parse(text = src)[[1]]
  target <- function(x) {
    if (!is.call(x)) return(x)
    nm <- if (is.name(x[[1]])) as.character(x[[1]]) else ""
    if (nm %in% c("frm", "frm_multiple")) {
      x[[arg]] <- str2lang(value)
      return(x)
    }
    for (i in seq_along(x)) if (is.call(x[[i]])) x[[i]] <- target(x[[i]])
    x
  }
  paste(deparse(target(e)), collapse = "\n")
}

AUTO_RETRY <- list(
  # brms defaults to gaussian(); frm() has no default family.
  "default-family" = function(src, msg) {
    if (!grepl("No family specified", msg)) return(NULL)
    add_arg(src, "family", "gaussian()")
  },
  # brms's update() takes formula.; frmtmb's takes formula.
  "update-formula-dot" = function(src, msg) {
    if (!grepl("unused argument \\(formula\\.", msg)) return(NULL)
    sub("formula\\.\\s*=", "formula =", src)
  },
  # brms's update() takes newdata; frmtmb's splices arguments straight
  # into the stored frm() call, where the argument is data.
  "update-newdata" = function(src, msg) {
    if (!grepl("unused argument \\(newdata", msg)) return(NULL)
    sub("newdata\\s*=", "data =", src)
  },
  # A family passed as a bare constructor name.
  "family-call" = function(src, msg) {
    if (!grepl("Cannot interpret `family` of class frmtmb_family", msg)) return(NULL)
    sub("family = ([a-zA-Z_.][a-zA-Z_.0-9]*)([,)])", "family = \\1()\\2", src)
  }
)

# Keyed by expression id. Value is the replacement source.
PATCH <- list(
  # brms lets one nlpar formula name several parameters; frmtmb needs one
  # formula per parameter.
  "brms_nonlinear.2.2" = 'fit1 <- frm(bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1, nl = TRUE),
    data = dat1, family = gaussian(), start = list(beta = c(1, 0.5)))',
  # Nonlinear fits need a starting region; in brms the priors supplied it.
  "brms_nonlinear.9.1" = 'fit_loss <- frm(bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
    ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE), data = loss,
    family = gaussian(), start = list(beta = c(5000, 1, 45)))',
  "brms_multilevel.15.3" = 'fit_loss1 <- frm(formula = nlform, data = loss,
    family = gaussian(), start = list(beta = c(5000, 1, 45)))',
  "brms_multilevel.19.2" = 'fit_loss2 <- update(fit_loss1, formula = nlform2,
    start = list(beta = c(5000, 1, 45)))',
  # update() takes a complete formula, not brms's one-sided delta.
  "brms_overview.7.1" =
    'fit2 <- update(fit1, formula = time | cens(censored) ~ age * sex + disease + (1 | patient))',
  "brms_phylogenetics.11.1" =
    'model_repeat2 <- update(model_repeat1,
       formula = phen ~ spec_mean_cf + within_spec_cf +
         (1 | gr(phylo, cov = A)) + (1 | species),
       data = data_repeat)',
  # lf() has no frmtmb equivalent; the dpar formula goes inside bf().
  "brms_multivariate.9.1" = 'bf_tarsus <- bf(tarsus ~ sex + (1|p|fosternest) + (1|q|dam),
    sigma ~ 0 + sex) + skew_normal()',
  # brms tests a one-sided hypothesis directly; frmtmb tests a contrast
  # against zero and the user reads the sign and the interval.
  "brms_distreg.5.1" = 'hyp <- "exp(sigma_Intercept + sigma_grouptreat) - exp(sigma_Intercept)"',
  "brms_overview.6.1" = 'hypothesis(fit1, "Intercept - age", class = "sd", group = "patient")'
)
