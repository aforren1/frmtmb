# Lane hollow: the bins 1 and 2 totals before and after this lane's
# reclassifications, read from the tracked ledgers with the rules of
# dev/brmsport-bin2-ledger.R (its ct_kind map and its totals). No
# recorded run exists in any tree, so the ledger script cannot be rerun
# here; the consolidation regenerates it.
#   Rscript dev/hollow-totals.R > dev/hollow-out/totals.txt 2>&1
rd <- function(f) utils::read.delim(f, quote = "", colClasses = "character",
                                    na.strings = NULL)
ct_kind <- c(stan = "decision", mcmc = "decision", "no-draws" = "decision",
             "brms-internal" = "decision", "no-density" = "decision",
             "no-quantile" = "decision", "brms-backends" = "decision",
             absent = "gap", "no-cdf" = "gap", "not-ported" = "cost",
             fixture = "artifact", upstream = "artifact",
             reference = "artifact", adapter = "artifact")
b1 <- rd("dev/brmsport-ledger.tsv")
b1$n <- 1L
b2 <- rd("dev/brmsport-bin2-ledger.tsv")
b2$n <- as.integer(b2$n)
stopifnot(nrow(b1) == 494L, sum(b2$n) == 329L)

totals <- function(b1, b2, label) {
  all <- rbind(b1[c("id", "outcome", "class", "n")],
               b2[c("id", "outcome", "class", "n")])
  kind <- ifelse(all$outcome == "cannot transfer",
                 ct_kind[all$class], "")
  stopifnot(!anyNA(kind))
  np <- sum(all$n[all$outcome == "pass"])
  nd <- sum(all$n[kind == "decision"])
  cat(sprintf("%-34s bins 1 and 2 pass %d of %d; without decisions %d of %d",
              label, np, sum(all$n), np, sum(all$n) - nd), "\n")
  b1p <- sum(b1$outcome == "pass")
  b1d <- sum(b1$outcome == "divergence")
  b1c <- sum(b1$outcome == "cannot transfer")
  cat(sprintf("%-34s bin 1 pass/divergence/cannot transfer %d/%d/%d",
              "", b1p, b1d, b1c), "\n")
  cat(sprintf("%-34s divergence (all bins) %d; hollow class %d", "",
              sum(all$n[all$outcome == "divergence"]),
              sum(all$n[all$class == "hollow"])), "\n")
}
show <- function(b, ids) {
  print(b[b$id %in% ids, c("id", "outcome", "class")], row.names = FALSE)
}
totals(b1, b2, "base (tracked ledgers):")

# Row 3: the two custom-family rows become decisions
r3 <- c("log_lik:603", "posterior_predict:440")
show(b2, r3)
b2$class[b2$id %in% r3] <- "brms-internal"
totals(b1, b2, "after Row 3:")

# Rows 1 and 2, as this lane classifies them (set below once decided)
row1 <- Sys.getenv("HOLLOW_ROW1", "")
row2 <- Sys.getenv("HOLLOW_ROW2", "")
set_row <- function(b, id, spec) {
  if (!nzchar(spec)) return(b)
  p <- strsplit(spec, ":", fixed = TRUE)[[1]]
  b$outcome[b$id == id] <- p[1]
  b$class[b$id == id] <- if (length(p) > 1L) p[2] else ""
  b
}
show(b1, "brmsfit-methods:394")
show(b2, "brmsfit-methods:469")
b1 <- set_row(b1, "brmsfit-methods:394", row1)
b2 <- set_row(b2, "brmsfit-methods:469", row2)
totals(b1, b2, "after Rows 1 and 2 as well:")
