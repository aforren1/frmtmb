# The verdict moves this lane makes to the ported brms tier
# (dev/brmsport-ledger.tsv), applied idempotently so the files can be
# rebuilt from the base tree.
#
#   Rscript dev/adefects-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# WHY A SCRIPT AND NOT AN EDIT. dev/brmsport-punch1-verdicts.R set the
# precedent: the three verdict files are data, and a hand edit to a
# 200-line TSV cannot be re-read or re-applied. Every move below names
# the row, what it was, what it becomes and the run that shows it
# (dev/adefects-log/port-record.txt).
man_f <- "dev/brmsport-verdicts-manual.tsv"
fit_f <- "dev/brmsport-verdicts-manual-fit.tsv"
own_f <- "dev/brmsport-verdicts-own.tsv"

rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
man <- rd(man_f)
fit <- rd(fit_f)
own <- rd(own_f)

# Rows that now HOLD, so their manual verdict is stale and the ledger
# counts them as passes. The ledger builder stops on a stale verdict, so
# these are the five it named.
to_pass <- c(
  # nrow(me[[2]]) on fit1$data, which is the 40-row frame now
  "brmsfit-methods:179",
  # hypothesis(fit3, "b_Age x 0"): frmtmb's refusal now CONTAINS brms's
  # own sentence, "Every hypothesis must be of the form 'left (= OR <
  # OR >) right'", so the row passes on brms's pattern as written
  "brmsfit-methods:417",
  # model.frame(fit1) == fit1$data, which is now an identity
  "brmsfit-methods:567",
  # pp_check(fit1, newdata = fit1$data[1:10, ])
  "brmsfit-methods:675",
  # is.null(attr(fit1$data$patient, "contrasts")): a real column and a
  # real attribute read now, where it used to hold on a NULL from the
  # data2 partial match and was hand-marked hollow
  "brmsfit-methods:1035")
fit <- fit[!(fit$id %in% to_pass), , drop = FALSE]

# Rows that still do not hold, but whose reason named the `fit$data`
# partial match, which is gone. Each reason is replaced by what the run
# now reports, and the class by what that fault is.
recl <- list(
  list(id = "brmsfit-methods:345", verdict = "defect", class = "shape",
       reason = paste0(
         "dim(fitted(fit1, newdata = fit1$data[1:10, ])) is NULL: ",
         "fitted() returns a vector where brms returns a 10 x 4 ",
         "summary matrix. Rule 3, item 2.6f; :348 is the same ",
         "assertion without newdata. Was fit-data, and the fit$data ",
         "blocker is gone (dev/adefects-findings.md D4)")),
  list(id = "brmsfit-methods:350", verdict = "defect", class = "shape",
       reason = paste0(
         "dim(fitted(fit4, newdata = fit4$data[1, ])) is 1 x 4 where ",
         "brms gives 1 x 4 x 4: an ordinal fitted() returns the ",
         "category probabilities, not brms's draws summary per ",
         "category. Rule 3, item 2.6f; :348 is the same shape without ",
         "newdata. Was fit-data")),
  list(id = "brmsfit-methods:352", verdict = "defect", class = "shape",
       reason = paste0(
         "dim(fitted(fit4, newdata = fit4$data[1, ], scale = ...)) is ",
         "NULL where brms gives 1 x 4 x 3. Rule 3, item 2.6f. Was ",
         "fit-data")),
  list(id = "brmsfit-methods:682", verdict = "defect",
       class = "internal-error",
       reason = paste0(
         "the assignment above it fails: pp_check(group = ) hands ",
         "bayesplot the group NAME rather than the column, so every ",
         "ppc_*_grouped type dies on bayesplot's 'length(group) must ",
         "be equal to the number of observations'. Measured with and ",
         "without newdata, and on a plain y ~ x + (1 | g) fit as well ",
         "as on fixture 1 (dev/adefects-findings.md, found and not ",
         "fixed 1). Was fit-data")),
  list(id = "brmsfit-methods:764", verdict = "pending 2.6d",
       class = "shape",
       reason = paste0(
         "dim(predict(fit4, newdata = fit4$data[1, ])) is NULL where ",
         "brms gives 1 x 4: predict() returns a vector, which item ",
         "2.6d is about; :761, :762 and :767 are the same assertion ",
         "on other rows. Was fit-data")),
  list(id = "brmsfit-methods:924", verdict = "defect", class = "output",
       reason = paste0(
         "attr(up$data, 'data_name') is NULL: brms records the name of ",
         "the newdata argument on the frame it stores and frmtmb does ",
         "not. The frame itself is there now. Was fit-data")),
  list(id = "data-helpers:7", verdict = "defect", class = "output",
       reason = paste0(
         "expect_silent() fails on R's 'contrasts dropped from factor ",
         "Trt' warning, twice, from model.matrix() on the newdata ",
         "slice; brms is silent. The NUMBERS are unaffected: the five ",
         "predictions equal the first five of the full-data ones, max ",
         "absolute difference 0 (dev/adefects-log/evidence.txt). Was ",
         "fit-data")),
  list(id = "data-helpers:12", verdict = "defect",
       class = "different-error",
       reason = paste0(
         "brms refuses the integer codes of `fac`, a factor column of ",
         "the data that the model does not use; frmtmb refuses the ",
         "new `visit` level set on the line before, and never checks ",
         "an unused column. Both refuse the call, for different ",
         "faults. Was fit-data")))
for (r in recl) {
  i <- which(fit$id == r$id)
  stopifnot(length(i) == 1L)
  fit$verdict[i] <- r$verdict
  fit$class[i] <- r$class
  fit$reason[i] <- r$reason
}

# Rows that now make brms's OWN refusal in frmtmb's words (rule 1).
# They move out of the manual files and into the own-words file, where
# the harness asserts frmtmb's pattern on brms's call.
own_new <- data.frame(
  id = c("brm:106", "brm:108", "data-helpers:9"),
  pkg = "frmtmb",
  pattern = c("the time index must be ONE variable name",
              "grouping term must be variable names combined by",
              "New levels in grouping factor .visit.: 5"),
  note = c(paste0("brms: Cannot coerce 'x + y' to a single variable ",
                  "name. The grammar check runs before the cov = TRUE ",
                  "one, so the ported call, which omits cov, reaches ",
                  "it (dev/adefects-log/evidence.txt)"),
           paste0("brms: Illegal grouping term 'g1/g2'. The grammar ",
                  "check runs before the cov = TRUE one, so the ported ",
                  "call reaches it"),
           paste0("brms: Levels '5' of grouping factor 'visit' cannot ",
                  "be used. The same refusal of the same case in ",
                  "frmtmb's words")),
  stringsAsFactors = FALSE)
man <- man[!(man$id %in% own_new$id), , drop = FALSE]
fit <- fit[!(fit$id %in% own_new$id), , drop = FALSE]
own <- own[!(own$id %in% own_new$id), , drop = FALSE]
own <- rbind(own, own_new)
own <- own[order(own$id), , drop = FALSE]

wr(man, man_f)
wr(fit, fit_f)
wr(own, own_f)
cat("manual", nrow(man), "manual-fit", nrow(fit), "own", nrow(own), "\n")
