# The verdict moves of the 0.72.0 consolidation in the ported brms tier,
# bin 1 (dev/brmsport-ledger.tsv), applied idempotently: a row is
# recorded only once every lane's code is in one build.
#
#   FRMTMB_PORT_LIB=<release lib> sh dev/defects-record.sh
#   Rscript dev/rel072-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# A pass carries no verdict (dev/brmsport-ledger.R stops on a verdict
# for an assertion that held). Lane prilabel: priors:74 ("default_prior
# returns correct prior names for auxiliary parameters") holds, because
# default_prior() now lists the per-coefficient sd row (sd, g,
# Intercept) that its divergence named (dev/prilabel-findings.md).
man_f <- "dev/brmsport-verdicts-manual.tsv"
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
wr <- function(x, f) utils::write.table(x, f, sep = "\t", quote = FALSE,
                                        row.names = FALSE)
man <- rd(man_f)
gone <- man$id %in% "priors:74" & man$pkg == "frmtmb"
cat("verdicts removed (now hold):", man$id[gone], "\n")
if (any(gone)) wr(man[!gone, , drop = FALSE], man_f)
