# The verdict moves of the 0.69.0 consolidation to the ported brms tier
# (dev/brmsport-ledger.tsv), applied idempotently. A row is recorded
# only once every lane's code is in one build, so the moves are made
# here, against the record of the merged tree.
#
#   FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r7 sh dev/defects-record.sh
#   Rscript dev/rel069-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# dev/round-20261007.md has the before and after counts;
# dev/rel069-log/ledger-stale.txt is the builder's list of the rows that
# hold, and dev/rel069-log/msgdiff.txt the rows whose message moved.
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

# Rows that now HOLD on brms's own assertion: lane optima's mo() fit
# reaches the maximum, so fixture 1 converges and fitted(fit1, dpar =
# 'sigma') has finite error columns (its review, claim 6)
to_pass <- c("brmsfit-methods:326")
# own-words rows whose brms assertion now holds as brms wrote it: lane
# surface's pp_mixture() of draws of a non-mixture model stops in brms's
# words (its review, m6)
own_pass <- c("brmsfit-methods:719")

stopifnot(sum(c(man$id, fit$id) %in% to_pass) <= length(to_pass))
man <- man[!(man$id %in% to_pass), , drop = FALSE]
fit <- fit[!(fit$id %in% to_pass), , drop = FALSE]
own <- own[!(own$id %in% own_pass), , drop = FALSE]
wr(man, man_f)
wr(fit, fit_f)
wr(own, own_f)
cat("removed", length(to_pass), "manual ids and", length(own_pass),
    "own-words ids\n")
