# The verdict moves of the 0.66.0 consolidation to the ported brms tier
# (dev/brmsport-ledger.tsv), applied idempotently after
# dev/defects-verdicts.R, which the release tree applied first. The six
# lanes of 2026-09-29 each moved rows with their own code; a row holds
# only once every lane's code is in one build, so the moves are made
# here, against the record of the merged tree.
#
#   FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r4 sh dev/defects-record.sh
#   Rscript dev/rel066-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# dev/round-20260929b.md has the before and after counts.
man_f <- "dev/brmsport-verdicts-manual.tsv"
fit_f <- "dev/brmsport-verdicts-manual-fit.tsv"

rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
man <- rd(man_f)
fit <- rd(fit_f)

# Rows that now HOLD on brms's own assertion: the ledger builder named
# each as a stale verdict on the merged build (dev/rel066-log/).
to_pass <- c(
  # formula2: dpar equations and cmc
  paste0("brmsformula:", c(30, 32, 34, 36, 38)), "standata:744",
  # aterms2: cat(), subset(), index(), mi(idx = ), rate()
  paste0("standata:", c(171, 266, 267, 268, 269, 620, 634, 641, 649,
                        1047)),
  # fams2: zero_inflated_beta_binomial(), hurdle_cumulative(), xbeta()
  paste0("families:", c(52:76, 81, 113)),
  # postfit2: conditional_smooths(), conditional_effects() options,
  # make_conditions(), update_adterms(), posterior_average()
  paste0("brmsfit-methods:", c(157, 170, 213, 215, 269, 275, 277, 619,
                               620)),
  "brmsfit-helpers:94", "brmsfit-helpers:95",
  paste0("brmsformula:", c(44, 48, 52, 56)),
  # formula2's family list with defects' default_prior() that no longer
  # reads the response: neither lane alone made these hold
  paste0("priors:", 91:94))

# Rows that still do not hold, whose recorded reason named a gap a lane
# of this round closed. Each reason is what the merged build reports.
recl <- list(
  list(id = "brmsfit-methods:685", verdict = "divergence",
       class = "no-draws",
       reason = paste0(
         "pp_check(type = 'loo_pit_qq') on a maximum-likelihood fit is ",
         "refused by name: a loo type weights posterior draws, and an ML ",
         "fit has none. frmtmb.sample's draws method builds the PSIS ",
         "weights as brms does (lane sampfix)")),
  list(id = "standata:928", verdict = "divergence", class = "policy",
       reason = paste0(
         "y ~ . expands against the data as brms's does (lane formula2); ",
         "brms's data has x1 == x2, and frmtmb drops the aliased x2 as ",
         "lm() does ('rank deficient; dropping column(s): x2'), where ",
         "brms keeps a column that only its prior places")),
  list(id = "brmsfit-methods:953", verdict = "divergence", class = "class",
       reason = paste0(
         "update(fit2, bf(. ~ ., a + b ~ 1, nl = TRUE)) refits now, ",
         "keeping the nonlinear body (lane formula2); the row asserts ",
         "is(up, 'brmsfit'), and frmtmb objects must NOT carry brms's ",
         "class names (user decision, 2026-09-17, rule 2), as ",
         "brmsfit-methods:957")),
  list(id = "brmsfit-methods:959", verdict = "divergence", class = "class",
       reason = paste0(
         "update(fit3, bf(~ ., family = acat())) refits now with family ",
         "acat and the old formula (lane formula2); the row asserts ",
         "is(up, 'brmsfit'), which rule 2 of 2026-09-17 forbids, as ",
         "brmsfit-methods:957")),
  list(id = "standata:738", verdict = "cannot transfer", class = "absent",
       reason = paste0(
         "lf(cmc = ) works now (lane formula2), but frmtmb's cumulative() ",
         "has no disc, so lf(disc ~ ...) is refused as a parameter the ",
         "family lacks; X_disc and Z_1_disc_1 are Stan data")),
  list(id = "standata:739", verdict = "cannot transfer", class = "absent",
       reason = paste0(
         "lf(cmc = ) works now, but frmtmb's cumulative() has no disc ",
         "(standata:738); Z_1_disc_1 is Stan data")),
  list(id = "standata:740", verdict = "cannot transfer", class = "absent",
       reason = paste0(
         "lf(cmc = ) works now, but frmtmb's cumulative() has no disc ",
         "(standata:738)")),
  list(id = "standata:745", verdict = "cannot transfer", class = "stan",
       reason = paste0(
         "bf(cmc = ) works now and frmtmb's Z equals brms's ",
         "(dev/formula2-findings.md 2.4); the row reads Z_1_1, which the ",
         "harness's standata view does not carry")),
  list(id = "brmsfit-methods:162", verdict = "defect", class = "argument",
       reason = paste0(
         "conditional_effects() takes too_far and surface now (lane ",
         "postfit2); plot() of the result refuses brms's plot = FALSE by ",
         "name, since frmtmb draws with base graphics where brms returns ",
         "ggplot objects (as brmsfit-methods:154 and :391)")),
  list(id = "brmsfit-methods:164", verdict = "defect", class = "argument",
       reason = paste0(
         "conditional_effects() takes too_far now (lane postfit2); plot() ",
         "refuses brms's stype and plot arguments by name (as ",
         "brmsfit-methods:162)")),
  list(id = "brmsfit-methods:167", verdict = "cannot transfer",
       class = "no-draws",
       reason = paste0(
         "conditional_effects() takes spaghetti now (lane postfit2); on a ",
         "maximum-likelihood fit it needs band = 'boot', one curve per ",
         "refit, and this call's default Wald band is refused by name: ",
         "the fit has only its estimate")),
  list(id = "brmsfit-methods:169", verdict = "cannot transfer",
       class = "no-draws",
       reason = paste0(
         "reads the spaghetti object of brmsfit-methods:167, refused on ",
         "an ML fit under the Wald band; plot() also refuses brms's ",
         "plot = FALSE by name")),
  list(id = "brmsfit-methods:273", verdict = "cannot transfer",
       class = "no-draws",
       reason = paste0(
         "conditional_smooths() exists now (lane postfit2); on a ",
         "maximum-likelihood fit it refuses ndraws and spaghetti by name, ",
         "since the fit has no draws")),
  list(id = "families:107", verdict = "cannot transfer",
       class = "brms-internal",
       reason = paste0(
         "xbeta() exists now (lane fams2); the row reads brms's family ",
         "field $closed, which no frmtmb family carries")),
  list(id = "families:108", verdict = "cannot transfer",
       class = "brms-internal",
       reason = paste0(
         "xbeta() exists now (lane fams2); the row reads brms's family ",
         "field $ybounds, which no frmtmb family carries")),
  list(id = "families:109", verdict = "cannot transfer",
       class = "brms-internal",
       reason = paste0(
         "xbeta() exists now (lane fams2); the row reads brms's family ",
         "field $type as 'real', where frmtmb's type is 'continuous'")),
  list(id = "families:117", verdict = "cannot transfer", class = "stan",
       reason = paste0(
         "xbeta() exists now (lane fams2); the row prints $prior, brms's ",
         "Stan prior string, which no frmtmb family carries")),
  list(id = "families:118", verdict = "cannot transfer", class = "stan",
       reason = paste0(
         "xbeta() exists now (lane fams2); the row prints $prior, brms's ",
         "Stan prior string (families:117)")),
  list(id = "families:122", verdict = "cannot transfer", class = "stan",
       reason = paste0(
         "xbeta() exists now (lane fams2); the row reads $include, the ",
         "Stan file brms's family includes")),
  list(id = "priors:14", verdict = "cannot transfer", class = "absent",
       reason = paste0(
         "sratio() takes no threshold = 'equidistant' (cse() exists now, ",
         "lane fams2); equidistant thresholds are not built ",
         "(dev/fams2-findings.md, 'Decided not to do')")),
  list(id = "brmsfit-methods:627", verdict = "cannot transfer",
       class = "brms-internal",
       reason = paste0(
         "posterior_average() exists now (lane postfit2); the setup line ",
         "reaches brms:::SW and is not run")),
  list(id = "brmsfit-methods:628", verdict = "cannot transfer",
       class = "brms-internal",
       reason = paste0(
         "posterior_average() exists now (lane postfit2); reads the ",
         "object of brmsfit-methods:627, whose setup is not run"))
)

man <- man[!(man$id %in% to_pass), , drop = FALSE]
fit <- fit[!(fit$id %in% to_pass), , drop = FALSE]
for (r in recl) {
  hit <- 0L
  for (nm in c("man", "fit")) {
    x <- get(nm)
    i <- which(x$id == r$id)
    if (!length(i)) next
    x$verdict[i] <- r$verdict
    x$class[i] <- r$class
    x$reason[i] <- r$reason
    assign(nm, x)
    hit <- hit + length(i)
  }
  if (!hit) stop("no verdict row for ", r$id)
}
wr(man, man_f)
wr(fit, fit_f)
cat("removed", length(to_pass), "manual ids; reclassified", length(recl),
    "\n")
