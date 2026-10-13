# The verdict moves of the 0.70.0 consolidation to the ported brms tier
# (dev/brmsport-ledger.tsv), applied idempotently. A row is recorded
# only once every lane's code is in one build, so the moves are made
# here, against the record of the merged tree.
#
#   FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r8 sh dev/defects-record.sh
#   Rscript dev/rel070-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# dev/round-20261013.md has the before and after counts;
# dev/rel070-log/ledger-stale.txt is the builder's stop on the row, and
# dev/rel070-log/msgdiff.txt the rows whose message moved.
fit_f <- "dev/brmsport-verdicts-manual-fit.tsv"
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
fit <- rd(fit_f)

# brmsfit-methods:326 no longer holds: fixture 1 stops at "false
# convergence (8)" again, where its sigma coefficients' sdreport()
# variances are negative. 0.69.0's pass was the luck of the stop point
# (lane optrap's review, item 6: from 40 jittered starts the row holds
# from 34 of 41 on rellib-r7 and fails from 7; the objective is
# unbounded within reach), and lane optrap's unit for this fixture's
# student response (8) moves the default stop. The verdict is 0.68.1's
# (4f5ea39f), with the reason brought up to date.
row <- data.frame(
  id = "brmsfit-methods:326", pkg = "frmtmb", verdict = "cannot transfer",
  class = "fixture",
  reason = paste0(
    "fitted(fit1, dpar = 'sigma') is brms's four-column summary, and ",
    "FIXTURE 1 DOES NOT CONVERGE (nlminb false convergence (8), an sd at ",
    "-723 to -1253 on the log scale, an objective unbounded within ",
    "reach), so its sigma coefficients' variances are negative, the ",
    "Est.Error and the two Q columns are NaN and all(fi > 0) is NA ",
    "rather than TRUE. The Estimate column is positive, which is what ",
    "the assertion is about; brms's fixture converged. It held at 0.69.0 ",
    "by the luck of the stop point (34 of 41 jittered starts; lane ",
    "optrap's review, item 6). A converged stand-in would transfer it ",
    "(dev/brmsport-findings.md section 3)"),
  stringsAsFactors = FALSE)
stopifnot(identical(names(row), names(fit)))
fit <- fit[fit$id != row$id, , drop = FALSE]
# in the file's order: by topic, then by line number
topic <- sub(":.*$", "", fit$id)
line <- suppressWarnings(as.integer(sub("^.*:", "", fit$id)))
at <- which(topic == "brmsfit-methods" & line > 326L)[1]
if (is.na(at)) at <- nrow(fit) + 1L
fit <- rbind(fit[seq_len(at - 1L), , drop = FALSE], row,
             fit[seq_len(nrow(fit) - at + 1L) + at - 1L, , drop = FALSE])
wr(fit, fit_f)
cat("verdicts: brmsfit-methods:326 is", row$verdict, "/", row$class,
    "at row", at, "\n")
