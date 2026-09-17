# Reviewer (brmsport): the ledger's totals after the review's
# reclassifications and the user's 2026-09-17 decisions. Reads
# dev/brmsport-ledger.tsv, changes nothing on disk, prints the tables.
#   Rscript dev/brmsport-rev-totals.R > dev/brmsport-log/rev-totals.txt
L <- utils::read.delim("dev/brmsport-ledger.tsv", quote = "",
                       stringsAsFactors = FALSE)
stopifnot(nrow(L) == 494L, !anyDuplicated(L$id))
re <- function(ids, outcome, class) {
  hit <- match(ids, L$id)
  if (anyNA(hit)) stop("unknown id: ", ids[is.na(hit)])
  L$outcome[hit] <<- outcome
  L$class[hit] <<- class
}
before <- table(factor(L$outcome))
# user rule 1: the same refusal in frmtmb's own words is a pass
re(c(paste0("brm:", c(75, 77, 79, 83, 85, 88, 90, 92)),
     paste0("brmsfit-methods:", c(203, 207, 330, 414, 719)),
     paste0("standata:", c(106, 109, 112, 121, 124, 127, 130, 179, 181,
                           184, 726, 1118)),
     "families:44", "brm:110"), "pass", "")
# wording rows that are NOT the same refusal
re("brm:81", "defect", "different-error")
re("families:102", "defect", "argument")
# hollow passes
re("brmsfit-methods:394", "divergence", "no-draws")
re("brmsfit-methods:698", "defect", "internal-error")
# user rule 3, and divergences whose citation does not say it
re(paste0("brmsfit-methods:", c(314, 317, 825, 887)), "defect", "shape")
re("brmsfit-methods:217", "defect", "output")
re(paste0("brmsfit-methods:", c(213, 215, 358, 359, 840, 841)), "defect",
   "refuses-accepted")
re(c("priors:55", "priors:59"), "defect", "naming")
re("brmsfit-methods:107", "cannot transfer", "absent")
# user rule 2: brms's class names are a recorded divergence
re(c("brmsfit-methods:373", "brmsfit-methods:270", "priors:27"),
   "divergence", "class")
# a feature frmtmb has under another spelling
re(paste0("brmsfit-methods:", c(995, 269, 275, 277)), "defect", "spelling")

lv <- c("pass", "defect", "divergence", "pending 2.6d", "cannot transfer")
cat("before:\n"); print(before[lv])
cat("after:\n"); print(table(factor(L$outcome, lv)))
cat("\nafter, outcome by class:\n")
print(as.data.frame(table(outcome = L$outcome, class = L$class),
                    stringsAsFactors = FALSE) |>
        subset(Freq > 0) |> (\(d) d[order(d$outcome, d$class), ])(),
      row.names = FALSE)
s <- L[L$tier %in% c("sample", "both"), ]
cat("\nfrmtmb.sample half (sample and both tiers):\n")
print(table(factor(s$outcome, lv)))
cat("sample-half rows:", nrow(s), "\n")
