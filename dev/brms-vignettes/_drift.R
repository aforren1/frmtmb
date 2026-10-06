# How the recorded edge labels have drifted from what the build does now.
#
#   BV_OUT=<dir> Rscript dev/brms-vignettes/_drift.R
#
# Why: each bv() call carries the label its author gave it at 0.42.0,
# and the label is part of the script, so a re-run cannot change it.
# What a re-run CAN show is a call whose outcome no longer fits its
# label. Four drifts are possible, and each one is a finding:
#
#   FAILS-NOW   labeled CLEAN, SPELLING or BEHAVIOR (all of which ran
#               when labeled) and errors today: a regression, or a
#               translation that a deliberate change made wrong.
#   RUNS-NOW    labeled MISSING or REFUSAL and runs today: the gap may
#               have closed. Read the logged value before believing it.
#   STILL-OPEN  labeled MISSING or REFUSAL and still errors.
#   AS-LABELED  everything else.
out <- Sys.getenv("BV_OUT", unset = tempdir())
# Only the vignette records: this script writes drift.csv into the same
# directory, and reading it back made a second run fail in rbind().
d <- do.call(rbind, lapply(list.files(out, "^brms.*[.]csv$", full.names = TRUE),
                           utils::read.csv, stringsAsFactors = FALSE))
d$edge[is.na(d$edge) | d$edge == ""] <- "CLEAN"
d <- d[d$kind != "data", ]
d$drift <- ifelse(
  d$edge %in% c("CLEAN", "SPELLING", "BEHAVIOR") & !d$ok, "FAILS-NOW",
  ifelse(d$edge %in% c("MISSING", "REFUSAL") & d$ok, "RUNS-NOW",
         ifelse(d$edge %in% c("MISSING", "REFUSAL"), "STILL-OPEN",
                "AS-LABELED")))
cat("## drift counts (rows = drift, columns = label)\n\n")
print(table(d$drift, d$edge))
cat("\n## drift by vignette\n\n")
print(table(d$vignette, d$drift))
for (k in c("FAILS-NOW", "RUNS-NOW", "STILL-OPEN")) {
  e <- d[d$drift == k, ]
  cat("\n## ", k, ": ", nrow(e), "\n\n", sep = "")
  for (i in seq_len(nrow(e))) {
    cat("- [", e$edge[i], "] ", e$vignette[i], " / ", e$label[i], "\n",
        sep = "")
    if (nzchar(e$why[i])) cat("    was: ", substr(e$why[i], 1, 200), "\n",
                              sep = "")
    if (!e$ok[i]) cat("    now: ", substr(gsub("\n", " ", e$msg[i]), 1, 300),
                      "\n", sep = "")
  }
}
utils::write.csv(d, file.path(out, "drift.csv"), row.names = FALSE)
