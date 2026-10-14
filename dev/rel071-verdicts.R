# The verdict moves of the 0.71.0 consolidation to the ported brms tier,
# bin 1 (dev/brmsport-ledger.tsv), applied idempotently. A row is
# recorded only once every lane's code is in one build, so the moves
# are made here, against the record of the merged tree.
#
#   FRMTMB_PORT_LIB=<release lib> sh dev/defects-record.sh
#   Rscript dev/rel071-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# A pass carries no verdict (dev/brmsport-ledger.R stops on a verdict
# for an assertion that held), so a row that now holds leaves the
# manual table. Lane postfix: brmsfit-methods:314 and :317 hold, because
# fitted() accepts brms's ndraws and changes nothing (a maximum-
# likelihood fit has no draws to average); their assertions are
# `dim(fi) == c(100, 4)`, which holds genuinely, not in a weaker form.
# brmsfit-methods:955 stays a defect (start values; the default-start
# error now names the model's parameters and the two ways to give
# starts, but offers none).
fit_f <- "dev/brmsport-verdicts-manual-fit.tsv"
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
fit <- rd(fit_f)
now_hold <- c("brmsfit-methods:314", "brmsfit-methods:317")
gone <- fit$id %in% now_hold & fit$pkg == "frmtmb"
cat("verdicts removed (now hold):", fit$id[gone], "\n")
fit <- fit[!gone, , drop = FALSE]
stopifnot("brmsfit-methods:955" %in% fit$id)
wr(fit, fit_f)
cat("brmsfit-methods:955 kept as",
    fit$verdict[fit$id == "brmsfit-methods:955"], "/",
    fit$class[fit$id == "brmsfit-methods:955"], "\n")
