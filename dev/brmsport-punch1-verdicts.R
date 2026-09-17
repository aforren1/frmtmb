# Punch round 1 (dev/reviews/20260917-brmsport.md): apply the review's
# reclassifications and the user's three rules of 2026-09-17 to the
# manual verdict files. Idempotent: every change is a set, not a toggle.
#
#   Rscript dev/brmsport-punch1-verdicts.R
files <- c("dev/brmsport-verdicts-manual.tsv",
           "dev/brmsport-verdicts-manual-fit.tsv")
rd <- function(f) {
  utils::read.delim(f, quote = "", colClasses = "character",
                    na.strings = NULL)
}
m <- lapply(files, rd)
names(m) <- files
own <- rd("dev/brmsport-verdicts-own.tsv")

set <- function(ids, verdict, class, reason = NULL, pkg = NULL,
                file = files[2]) {
  for (id in ids) {
    hit <- FALSE
    for (f in files) {
      i <- which(m[[f]]$id == id)
      if (length(i)) {
        m[[f]]$verdict[i] <<- verdict
        m[[f]]$class[i] <<- class
        if (!is.null(reason)) m[[f]]$reason[i] <<- reason
        hit <- TRUE
      }
    }
    if (!hit) {
      stopifnot(!is.null(reason), !is.null(pkg))
      m[[file]] <<- rbind(m[[file]], data.frame(
        id = id, pkg = pkg, verdict = verdict, class = class,
        reason = reason))
    }
  }
}
drop <- function(ids) {
  for (f in files) m[[f]] <<- m[[f]][!m[[f]]$id %in% ids, ]
}

rule2 <- paste("brms's class name; frmtmb objects must NOT carry brms's",
               "class names (user decision, 2026-09-17, rule 2)")
rule3 <- "(user decision, 2026-09-17, rule 3; item 2.6f)"

# rule 1: the 27 rows in dev/brmsport-verdicts-own.tsv are passes in
# frmtmb's own words and leave the manual files
drop(own$id)

set("brm:81", "defect", "different-error", paste(
  "not the same refusal: frmtmb evaluates the data first and stops on",
  "object 'sei' not found, where brms refuses se() for weibull from the",
  "formula alone; with sei present frmtmb refuses se() for weibull",
  "(dev/brmsport-probe3.R)"))
set("families:102", "defect", "argument", paste(
  "mixture() refuses the argument NAME order, where brms refuses its",
  "value 'x' (famlink ledger row 77, not fixed)"))

# MAJOR 4: hollow passes, now rejected by the harness's own rules
set("brmsfit-methods:394", "divergence", "no-draws", paste(
  "Evid.Ratio is NA on a fit (D4, dev/brmsnames-findings.md 'The",
  "hypothesis() object'); brms's assertion held only as",
  "is.numeric(NA_real_), which the harness now rejects as a type check",
  "on an all-NA value"), pkg = "frmtmb")
set("brmsfit-methods:698", "defect", "internal-error", paste(
  "pp_check(fit2, 'violin_grouped') dies with R's 'argument group is",
  "missing, with no default'; brms refuses by design ('Argument group is",
  "required'); brms's regex 'group' matched the internal error, which the",
  "harness now rejects"), pkg = "frmtmb")

# rule 3 and review MAJOR 6/7
set(c("brmsfit-methods:314", "brmsfit-methods:317"), "defect", "shape",
    paste("fitted() refuses ndraws, and past the refusal it returns a",
          "vector where brms returns the 4-column summary", rule3))
set("brmsfit-methods:825", "defect", "shape",
    paste("residuals() refuses probs, and past the refusal it returns a",
          "vector where brms returns the summary", rule3))
set("brmsfit-methods:887", "defect", "shape",
    paste("summary() refuses priors = TRUE, and past the refusal",
          "summary(fit) has no $fixed or $random", rule3))
set(c("brmsfit-methods:358", "brmsfit-methods:359"), "defect",
    "refuses-accepted",
    paste("fitted() exists and refuses the multivariate fit6 as 'not",
          "supported yet'; brms answers", rule3))
set(c("brmsfit-methods:840", "brmsfit-methods:841"), "defect",
    "refuses-accepted",
    paste("residuals() exists and refuses the multivariate fit6 as 'not",
          "supported yet'; brms answers", rule3))
for (f in files) {
  i <- m[[f]]$class %in% c("shape", "output") &
    !grepl("rule 3", m[[f]]$reason) & m[[f]]$verdict == "defect"
  m[[f]]$reason[i] <- paste(m[[f]]$reason[i], rule3)
}
set("brmsfit-methods:217", "defect", "output", paste(
  "no 'Predictions are treated as continuous variables' warning on the",
  "ordinal fit4; the brmsnames note cited before is about draws, and the",
  "fit method was never audited (brmsnames 'Found and NOT fixed' 11)"))
set(c("brmsfit-methods:213", "brmsfit-methods:215"), "defect",
    "refuses-accepted", paste(
      "conditional_effects() refuses its default Wald band on the",
      "nonlinear fit2 and names band = 'boot'; brms answers. No document",
      "records refusing the default band as a decision"))
set(c("priors:55", "priors:59"), "defect", "naming", paste(
  "the prior table names a smooth's unpenalized column s(x).fx1 where",
  "brms writes sx_1, and variables() already writes bs_sx_1; priorform",
  "P11/P12 deferred to brmsnames and brmsnames deferred back, so no",
  "document decides it"), file = files[1])
set("brmsfit-methods:107", "cannot transfer", "absent", paste(
  "as.mcmc(inc_warmup = TRUE): frm_sample() stores no warmup, a recorded",
  "gap (dev/brms-api-diff.md (c), 'Blocked, not small'), not a decision"))

# rule 2
for (f in files) {
  i <- m[[f]]$class == "class"
  m[[f]]$reason[i] <- sub("[(]D9[^)]*[)]|[(]divergence list D9[^)]*[)]",
                          "", m[[f]]$reason[i])
  m[[f]]$reason[i] <- paste(trimws(m[[f]]$reason[i]), rule2)
}
set("brmsfit-methods:373", "divergence", "class",
    paste("is.brmsformula(formula(fit1)):", rule2))
set("brmsfit-methods:270", "divergence", "class",
    paste("is(ms, 'brms_conditional_effects'):", rule2))
set("priors:27", "divergence", "class",
    paste("expect_is(bprior, 'brmsprior'): the object is a",
          "frmtmb_priorlist;", rule2))

# MAJOR 8: spelling defects
set("brmsfit-methods:995", "defect", "spelling", paste(
  "parnames() is brms's live spelling of variables(); frmtmb.sample",
  "defines parnames() and refuses it, core has none",
  "(dev/brms-suite-audit.md section 7 contract 14)"))
set(c("brmsfit-methods:269", "brmsfit-methods:275",
      "brmsfit-methods:277"), "defect", "spelling", paste(
  "conditional_smooths() is absent, though R/conditional-effects.R says",
  "conditional_effects() 'also covers what brms calls",
  "conditional_smooths()'"))

# MAJOR 1
set("priors:74", "defect", "refuses-accepted", paste(
  "default_prior() and get_prior() validate the response against the",
  "family and refuse brms's rnorm response under Beta(); brms returns 7",
  "rows. PRE-EXISTING: refused at every commit back to base 0.58.0",
  "(dev/brmsport-rev-regress.R). priorform's P13 showed a table only",
  "because dev/priorform-ledger.R replaced brms's y with runif(10)"),
  file = files[1])

# MAJOR 3: the fit-data mechanism
for (f in files) {
  i <- m[[f]]$class %in% c("fit-data", "hollow")
  m[[f]]$reason[i] <- gsub(
    "a frmtmb fit has no data slot: fit\\$data is an empty list|a frmtmb fit has no data slot|which a frmtmb fit does not have|fit1\\$data is an empty list",
    "fit$data partial-matches fit$data2 (the fit has no data element)",
    m[[f]]$reason[i])
  m[[f]]$reason[i] <- paste(
    m[[f]]$reason[i],
    "[fit$data is the $ partial match of fit$data2, dev/brmsport-rev-silent.R R8]")
}

# MINOR 10, settled in punch round 2 (recheck R5, the main session's
# decision): frmtmb's compatibility table records the nonlinear and the
# multivariate emmeans cases as refused for now (R/compat.R:1356 and
# :1612, enforced by emm_mu_linpred()), so they are a declared limitation:
# cannot transfer, absent. The mo() rows :21 and :24 are declared
# nowhere and stay defects.
set(paste0("emmeans:", c(27, 35, 38, 50)), "cannot transfer", "absent",
    paste("frmtmb's emmeans support is declared refused for a nonlinear mu",
          "('needs a linear mu predictor', R/compat.R:1612, enforced by",
          "emm_mu_linpred()); the reason is hidden by emmeans, a",
          "separate defect"))
set("emmeans:42", "cannot transfer", "absent",
    paste("frmtmb's emmeans support is declared univariate only",
          "(R/compat.R:1356); the reason is hidden by emmeans, a separate",
          "defect"))

# punch round 2 (recheck R3/R4): the all-NA type-check rule was removed
# because it also rejects an NA that is the correct answer, so row 394
# is marked hollow by hand; the partial-$ rule now rejects row 1035
# automatically, so it is an ordinary fit-data defect
set("brmsfit-methods:394", "divergence", "hollow", paste(
  "holds only as is.numeric(NA_real_): Evid.Ratio is NA on a fit (D4,",
  "dev/brmsnames-findings.md 'The hypothesis() object'), where brms",
  "reports a number. Marked by hand, because a rule cannot tell this NA",
  "from an NA that is the right answer"), pkg = "frmtmb")
set("brmsfit-methods:1035", "defect", "fit-data")

# silent rows, reasons corrected by the review (MAJOR 2, MINOR 1, 3)
set("brm:106", "defect", "silent", paste(
  "refused only because cov = TRUE is missing; with cov = TRUE an",
  "expression time term FITS with no refusal: ar(x + t, g) logLik",
  "-30.33421 and ar(t - 10 * x, g) -31.00618 against -30.88513 for",
  "ar(t, g), and ma(x + t, g) -29.84347 against -30.49819; brms refuses",
  "all three (dev/brmsport-defects.R S2)"))
set("brmsfit-methods:713", "defect", "silent", paste(
  "posterior_epred() ignores brms's point_estimate = 'median' and",
  "ndraws_point_estimate = 2 and returns all 25 draws, 25 x 40, where",
  "brms returns 2 x 40. Ignoring an UNKNOWN argument is brms parity:",
  "brms also returns 25 x 40 for not_an_argument = 2",
  "(dev/brmsport-defects.R S1)"))
set("brmsfit-methods:714", "defect", "silent", paste(
  "as brmsfit-methods:713: the point estimate brms repeats is never",
  "computed, so the rows differ"))
set("brmsfit-methods:417", "defect", "misparse", paste(
  "'b_Age x 0' is refused loudly as the unknown parameter b_b_Agex0,",
  "not as malformed. The silent half is elsewhere: a hypothesis with no",
  "relation, 'Trt1 + Age', returns exactly the row of 'Trt1 + Age = 0'",
  "where brms refuses (dev/brmsport-defects.R S5)"))
set("brmsfit-methods:837", "defect", "accepts-refused", paste(
  "residuals() on the ordinal fit4 returns the response code minus the",
  "expected score, y - sum(k P(Y = k)), identical at relative 0; brms",
  "refuses because predictive errors are not defined for ordinal",
  "models (dev/brmsport-defects.R S7)"))
set("brmsfit-methods:888", "defect", "shape", paste(
  "blocked by the priors refusal, and past it summary(fit) has no",
  "$fixed: its slots are frmtmb's (call, family, ..., coefficients,",
  "varcor, rescor, autocor, smooth_edf, extras, fixed_dpars)", rule3))

for (f in files) {
  write.table(m[[f]], f, sep = "\t", quote = FALSE, row.names = FALSE)
}
all <- do.call(rbind, m)
print(table(all$verdict))
