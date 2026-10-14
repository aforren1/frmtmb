# The verdicts of bin 2T, the 59 assertions that read a posterior, and
# the object each is judged on. Writes dev/brmsport-bin2-targets.tsv and
# dev/brmsport-bin2-verdicts-manual.tsv, which dev/brmsport-bin2-ledger.R
# reads. Written as a script so the reasons are reviewable R strings and
# not hand-typed tab-separated text.
#
#   Rscript dev/brmsport-bin2-verdicts-2t.R
#
# THE TARGET. A 2T assertion reads a posterior. frmtmb.sample's draws
# are frmtmb's posterior, so the 58 methods rows are judged on
# frm_sample() draws of the fixture (target frmtmb.sample), brms's text
# unchanged. The same text also runs on the maximum-likelihood fixture
# fit (frmtmb); that run is the restatement for ML, recorded beside the
# verdict and not counted. priors:130 is core-tier (default_prior()),
# judged on frmtmb. Constructions: dev/remeasure-bin2-probe1.R, log
# dev/remeasure-out/probe/bin2-probe1.txt.
m <- "brmsfit-methods"
ids <- function(...) paste0(m, ":", c(...))
targets <- rbind(
  data.frame(id = ids(119, 121, 123, 136, 138, 140, 142, 433, 434, 435,
                      442, 443, 446, 447, 448, 449, 452, 455, 457, 460,
                      467, 469, 471, 478, 479, 480, 487, 490, 497, 508,
                      511, 512, 519, 522, 523, 526, 533, 544, 547, 548,
                      555, 562, 563, 964, 966, 968, 969, 971, 972, 974,
                      1005, 1010, 1015, 1016, 1018, 1020, 1021, 1024),
             target = "frmtmb.sample"),
  data.frame(id = "priors:130", target = "frmtmb"))
stopifnot(nrow(targets) == 59L, !anyDuplicated(targets$id))

v <- list()
add <- function(id, pkg, verdict, class, reason) {
  v[[length(v) + 1L]] <<- data.frame(id = id, pkg = pkg, verdict = verdict,
                                     class = class, reason = reason)
}
S <- "frmtmb.sample"
F <- "frmtmb"

## On the draws: the verdicts that count --------------------------------
add(ids(121), S, "defect", "argument", paste(
  "bayes_R2() on draws refuses newdata and re_formula by name ('has no",
  "argument'); brms's bayes_R2() passes both to posterior_epred(), and",
  "frmtmb.sample's posterior_epred() takes newdata"))
add(ids(142), S, "cannot transfer", "fixture", paste(
  "coef(fit4 draws)$subject is 8 x 4 x 7: brms's 8 coefficients less",
  "disc_Intercept, because fixture 4 drops disc ~ 1 (frmtmb's ordinal",
  "families have no disc). The other 7 names are brms's"))
add(ids(434), S, "defect", "method", paste(
  "logLik() of draws has no method ('no applicable method'); brms's",
  "logLik() of a brmsfit is log_lik(), the draws x observations matrix,",
  "25 x 40 on brms's fixture 1"))
for (k in c(442, 443)) add(ids(k), S, "divergence", "policy", paste(
  "LOO() is brms's deprecated spelling; frmtmb.sample refuses it by name",
  "and names loo()"))
for (k in 446:449) add(ids(k), S, "defect", "argument", paste(
  "loo(fit1, fit1): frmtmb.sample's loo() takes one model and refuses",
  "several by name, naming loo_compare(); brms returns a loolist with",
  "the comparison. Not a decision: dev/brms-api-diff.md (c) lists",
  "compare = and model_names = on loo() and waic() as Small open work",
  "(fix round 1; owner lane postfix)"))
add(ids(469), S, "defect", "hollow", paste(
  "holds on a NULL newdata: draws carry no $data, so fit6$data is NULL,",
  "and loo() on draws swallows newdata in its dots anyway: loo(d6,",
  "newdata = 20 rows) returns all 40 pointwise rows",
  "(dev/remeasure-bin2-probe1.R). brms computes the elpd on newdata.",
  "Marked by hand"))
add(ids(471), S, "cannot transfer", "upstream", paste(
  "loo 2.10.1's loo_compare() returns a data frame whose first column is",
  "the model name, so [2, 1] is a string and the range check is FALSE;",
  "brms 2.23.0's own fixture 6 gives the same string",
  "(dev/remeasure-bin2-probe1.R)"))
for (k in 478:480) add(ids(k), S, "divergence", "policy", paste(
  "loo_subsample() on draws is refused by name: its point is a cheap",
  "approximation for data too large to evaluate, which a frmtmb model",
  "is not"))
absent <- list(c(487, 490, "loo_R2"), c(497, 508, 511, 512, "loo_epred"),
               c(519, 522, 523, 526, "loo_linpred"),
               c(533, 544, 547, 548, "loo_predict"),
               c(555, "loo_predictive_interval"),
               c(562, 563, "loo_model_weights"))
for (a in absent) {
  fn <- a[length(a)]
  for (k in a[-length(a)]) add(ids(k), S, "cannot transfer", "absent",
    paste0(fn, "() does not exist in frmtmb or frmtmb.sample ('could not ",
           "find function'). On draws it has a meaning (brms computes it ",
           "from the PSIS weights loo() already builds), so this is a ",
           "gap, not a decision"))
}
add(ids(1005), S, "divergence", "policy", paste(
  "WAIC() is brms's deprecated spelling; frmtmb.sample refuses it by",
  "name and names waic()"))
for (k in c(1015, 1016)) add(ids(k), S, "defect", "argument", paste(
  "waic(fit1, fit1): frmtmb.sample's waic() takes one model and refuses",
  "several by name, as loo() does (brmsfit-methods:446): open work in",
  "dev/brms-api-diff.md (c), owner lane postfix"))
add(ids(1021), S, "divergence", "policy", paste(
  "reads waic1, the result of WAIC(), refused by name",
  "(brmsfit-methods:1005); compare_ic() is brms's deprecated comparison",
  "and frmtmb has none"))
add("priors:130", F, "divergence", "policy", paste(
  "default_prior() describes frm(), which is flat in every slot (the",
  "policy of bin-1 priors:85); brms's student_t(3, -8, 2.5) is its",
  "sampling default for the intercept, shifted by the offset"))

## On the ML fit: the restatement, recorded and not counted -------------
nd <- function(k, why) add(ids(k), F, "divergence", "no-draws", why)
for (k in c(119, 121, 123)) nd(k, paste(
  "bayes_R2() is refused by name on a maximum-likelihood fit: it is",
  "computed per posterior draw"))
add(ids(136), F, "divergence", "policy", paste(
  "coef() of a fit is levels x coefficients, lme4's shape, without",
  "brms's Est.Error and quantile layer, the policy of bin-1",
  "brmsfit-methods:811; it also omits brms's sigma_Intercept and",
  "sigma_Trt1 columns (4 x 7 against 4 x 4 x 9)"))
nd(138, "coef(summary = FALSE) is refused by name on a fit (D3)")
add(ids(140), F, "divergence", "policy", paste(
  "coef() of a nonlinear fit is nested by nonlinear parameter, coef(fit)$a",
  "and $b, each levels x coefficients (the policy of brmsfit-methods:136),",
  "so coef(fit2)$patient is NULL; brms has one patient array",
  "10 x 4 x 4"))
add(ids(142), F, "divergence", "policy", paste(
  "coef() of a fit is levels x coefficients (brmsfit-methods:136),",
  "8 x 8 here; the fixture also drops disc ~ 1 (brmsfit-methods:142 on",
  "draws)"))
for (k in c(433, 435)) nd(k, paste(
  "log_lik() is refused by name on a fit: it is the pointwise likelihood",
  "at every posterior draw"))
nd(434, paste("ndraws() is refused by name on a fit; logLik() of a fit is",
              "stats' scalar log-likelihood"))
for (k in c(442, 443, 446:449, 452, 455, 457, 460, 467, 469, 471)) nd(k,
  paste("loo() and LOO() are refused by name on a maximum-likelihood fit,",
        "naming AIC() and loo(frm_sample(fit))"))
for (k in c(478:480, 487, 490, 497, 508, 511, 512, 519, 522, 523, 526,
            533, 544, 547, 548, 555, 562, 563)) nd(k, paste(
  "a leave-one-out quantity over posterior draws, which a",
  "maximum-likelihood fit has none of; the function is also absent",
  "from frmtmb core"))
nd(971, "VarCorr(summary = FALSE) is refused by name on a fit (D3)")
for (k in c(1005, 1010, 1015, 1016, 1018, 1020, 1021, 1024)) nd(k, paste(
  "waic(), WAIC() and add_criterion('waic') are refused by name on a",
  "maximum-likelihood fit: WAIC averages over posterior draws"))

V <- do.call(rbind, v)
# the verdict moves of the 0.71.0 consolidation
V <- source("dev/rel071-bin2-moves.R")$value(V, "2t")
stopifnot(!anyDuplicated(V[c("id", "pkg")]),
          all(V$id %in% targets$id), all(nzchar(V$reason)),
          !any(grepl("\t", V$reason)))
utils::write.table(targets, "dev/brmsport-bin2-targets.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
utils::write.table(V, "dev/brmsport-bin2-verdicts-manual.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("targets:", nrow(targets), " manual verdicts:", nrow(V), "(",
    sum(V$pkg == S), "on draws,", sum(V$pkg == F), "on the fit )\n")
