# The verdict moves of the 0.68.0 consolidation to the ported brms tier
# (dev/brmsport-ledger.tsv), applied idempotently. A row is recorded
# only once every lane's code is in one build, so the moves are made
# here, against the record of the merged tree.
#
#   FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r6 sh dev/defects-record.sh
#   Rscript dev/rel068-verdicts.R
#   Rscript dev/brmsport-ledger.R
#   Rscript dev/brmsport-gen.R
#
# dev/round-20261005.md has the before and after counts;
# dev/rel068-log/ledger-stale.txt is the builder's list of the 3 rows
# that hold, and dev/rel068-log/msgdiff.txt the rows whose message
# moved.
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

# Rows that now HOLD on brms's own assertion: lane gpby's gp(by = ) and
# the kriging covariance in emmeans() (its review: on the merged tree
# only these two flip)
to_pass <- c("brmsfit-methods:991", "emmeans:42")
# own-words rows whose brms assertion now holds as brms wrote it: cs()
# on a family that takes none is refused in brms's words since the
# consolidation's lift of cs() on cumulative() (small change a)
own_pass <- c("brm:116")

recl <- list(
  # the reason the fixes review checked (its final check, item 9;
  # dev/fixes-rev3-u955.R): from zero starts the ported update stops,
  # from intercepts at mean(count)/2 it reaches the ridge, and with
  # fit2's priors it converges
  list(id = "brmsfit-methods:955", verdict = "defect", class = "start-values",
       reason = paste0(
         "update(fit2, formula. = bf(count ~ a + b, nl = TRUE)) keeps ",
         "fit2's parameter formulas with brms's message, but the fixture ",
         "lacks fit2's priors; the update stops at its zero start ",
         "('NA/NaN gradient evaluation'); with the priors it converges ",
         "(7.1488, 11.5849, 1.6474, -0.7932; dev/fixes-rev3-u955.R). ",
         "brms's testmode does not fit, so brms never meets the start"))
)

man <- man[!(man$id %in% to_pass), , drop = FALSE]
fit <- fit[!(fit$id %in% to_pass), , drop = FALSE]
own <- own[!(own$id %in% own_pass), , drop = FALSE]
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
wr(own, own_f)
cat("removed", length(to_pass), "manual ids and", length(own_pass),
    "own-words ids; reclassified", length(recl), "\n")
