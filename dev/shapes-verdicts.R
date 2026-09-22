# Items 2.6d and 2.6f: move the ported tier's verdicts to what the code
# now does. Idempotent - every change is a set or a drop, not a toggle.
#
#   Rscript dev/shapes-verdicts.R
#   Rscript dev/brmsport-ledger.R      # recomputes and refuses a stale one
#   Rscript dev/brmsport-gen.R         # regenerates the tier
#
# A row is DROPPED when the assertion now holds: `dev/brmsport-ledger.R`
# takes a pass from the RUN and refuses a manual verdict that says
# otherwise, so a fixed defect is a deleted row and never a typed
# "pass".

files <- c("dev/brmsport-verdicts-manual.tsv",
           "dev/brmsport-verdicts-manual-fit.tsv")
rd <- function(f) {
  utils::read.delim(f, quote = "", colClasses = "character",
                    na.strings = NULL)
}
m <- lapply(files, rd)
names(m) <- files

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

# ---- fixed: the assertion now holds, so the verdict row goes ---------

# fitted() is brms's four-column summary, and the n x 4 x K array on an
# ordinal fit; it also takes nlpar and allow_new_levels now
drop(paste0("brmsfit-methods:", c(292, 293, 325, 328, 334, 337, 339,
                                  340, 342, 348, 355)))
# fixef() is brms's matrix, in brms's row order, and takes brms's pars
drop(paste0("brmsfit-methods:", c(364, 369)))
# ngrps() is brms's named list
drop(paste0("brmsfit-methods:", c(584, 585)))
# print(fit) is brms's layout, and summary() carries $fixed and $random
# and honors priors = TRUE
drop(paste0("brmsfit-methods:", c(784, 887, 888, 894, 896, 897)))
# residuals() is brms's four-column summary and takes probs
drop(paste0("brmsfit-methods:", c(825, 835)))
# vcov() is brms's 9 x 9 population-level block and takes correlation
drop(paste0("brmsfit-methods:", c(999, 1000)))
# predict() is brms's predictive summary (item 2.6d)
drop(paste0("brmsfit-methods:", c(726, 727, 729, 750, 753, 758, 761,
                                  762, 767)))
# frmtmb.sample: nsamples() and posterior_samples() answer as brms's do
drop(paste0("brmsfit-methods:", c(593, 594, 633, 634)))
# frmtmb.sample: point_estimate and ndraws_point_estimate are honored
drop(paste0("brmsfit-methods:", c(713, 714)))

# ---- reclassified: still not a pass, for a different reason ---------

set("brmsfit-methods:326", "cannot transfer", "fixture", paste(
  "fitted(fit1, dpar = 'sigma') is brms's four-column summary now, and",
  "FIXTURE 1 DOES NOT CONVERGE (nlminb code 1, NaN standard errors),",
  "so the Est.Error and the two Q columns are NA and all(fi > 0) is NA",
  "rather than TRUE. The Estimate column is positive, which is what the",
  "assertion is about; brms's fixture converged. A converged stand-in",
  "would transfer it (dev/brmsport-findings.md section 3)"),
    pkg = "frmtmb")

set("brmsfit-methods:891", "divergence", "no-draws", paste(
  "summary(fit)$fixed carries brms's four columns and then the Wald",
  "test this package reports, where brms writes Rhat, Bulk_ESS and",
  "Tail_ESS. Those three describe a SAMPLER and a maximum-likelihood",
  "fit has none, so the column set cannot match (item 2.6f)"))

set("brmsfit-methods:595", "defect", "argument", paste(
  "nsamples() answers now, as brms's does, but incl_warmup = TRUE is",
  "refused: frm_sample() discards the warmup rather than storing it, so",
  "there is nothing to count. The gap is dev/brms-api-diff.md (c),",
  "'Blocked, not small'"))

set("brmsfit-methods:635", "defect", "naming", paste(
  "posterior_samples(pars = '^b_') answers now and returns the right",
  "SET of columns; their ORDER is frmtmb's variables() order, which",
  "lists each predictor's coefficients together, where brms lists every",
  "intercept first. fixef(), vcov() and summary()$fixed take brms's",
  "order at item 2.6f; variables() keeps its own, and reordering it is",
  "a separate decision (dev/shapes-findings.md section 8)"))

# punch round 1: sample_new_levels = "gaussian" is ANSWERED now, so
# :775 runs and fails for a DIFFERENT, pre-existing reason (the
# newdata is fit5$data, which partial-matches fit5$data2 and has 2
# rows, the same silent defect :764 records); "old_levels" is still
# refused, so :772 stays a defect whose setup line dies
set("brmsfit-methods:772", "defect", "argument", paste(
  "predict() honors sample_new_levels = \"gaussian\", which is what it",
  "does, and refuses \"old_levels\" BY NAME: that value resamples the",
  "POSTERIOR draws of the levels the fit saw, and a",
  "maximum-likelihood fit has no such draws. The setup line dies, so",
  "this assertion reads a stale object"))

set("brmsfit-methods:775", "defect", "fit-data", pkg = "frmtmb",
    reason = paste(
  "predict(sample_new_levels = \"gaussian\") is ANSWERED now: an unseen",
  "level's effect is drawn from its block's estimated covariance. The",
  "assertion still fails, on the pre-existing fixture defect :764",
  "records: newdata is fit5$data[1:5, ], fit$data partial-matches",
  "fit$data2 (the fit has no data element), so newdata has 2 rows and",
  "the answer is 2 x 4 [fit$data is the $ partial match of fit$data2,",
  "dev/brmsport-rev-silent.R R8]"))

set(paste0("brmsfit-methods:", c(314, 317)), "defect", "argument", paste(
  "fitted() is brms's four-column summary now (item 2.6f) and still",
  "refuses ndraws BY NAME: a maximum-likelihood fit has no draws to",
  "thin. The setup line dies on that refusal, so the assertion reads a",
  "stale object"))

set("brmsfit-methods:747", "defect", "argument", paste(
  "predict(ndraws =) is answered now - the draws are SIMULATED, so",
  "ndraws sets how many - and the setup still dies earlier, on",
  "fit1$data, which partial-matches fit1$data2 (the fit has no data",
  "element)"))

for (f in files) {
  write.table(m[[f]], f, sep = "\t", quote = FALSE, row.names = FALSE)
}
all <- do.call(rbind, m)
cat("rows:", nrow(all), "\n")
print(table(all$verdict))
