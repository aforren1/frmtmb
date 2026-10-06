# Reviewer, re-check: the spell pass after the two patch deletions.
# Compares the reviewer's new spell pass with (a) the lane's new spell
# pass and (b) the reviewer's spell pass under the OLD patches.R, and
# rescores classes with classify() on the new pass.
#
#   Rscript dev/vigport-rev2-spellcmp.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
rd <- function(d) {
  out <- list()
  for (f in list.files(d, "[.]rds$", full.names = TRUE)) {
    out <- c(out, readRDS(f))
  }
  out
}
st <- function(x) vapply(x, function(r) r$status, "")
new <- rd(file.path(root, "vigport-rev2-out/spell/results-spell"))
lane <- rd(file.path(root, "vigport-port-out/r5/results-spell"))
old <- rd(file.path(root, "vigport-rev-out/seeded/results-spell"))
kind <- vapply(new, function(r) r$kind, "")
k <- names(new)
cat("rows", length(new), length(lane), length(old), "\n")
cat("new vs lane: status differs", sum(st(new)[k] != st(lane)[k]), "\n")
d <- k[st(new)[k] != st(old)[k]]
cat("new vs old patches: status differs", length(d), "\n")
for (id in d) cat(sprintf("  %-22s %-5s %s -> %s\n", id, kind[id],
                          st(old)[id], st(new)[id]))
for (kk in c("model", "post", "other")) {
  s <- k[kind == kk & st(new)[k] != "SETUP-SKIP"]
  cat(sprintf("%-5s spell OK old %d -> new %d of %d\n", kk,
              sum(st(old)[s] == "OK"), sum(st(new)[s] == "OK"), length(s)))
}
# classes: rebuild the merged table the way summarize.R does
HERE <- file.path(root, "brms-port")
source(file.path(HERE, "port-lib.R"))
M <- readRDS(file.path(root, "vigport-port-out/r5/results-merged.rds"))
M2 <- M
M2$status_spell <- st(new)[M2$id]
M2$msg_spell <- vapply(new[M2$id], function(r) r$msg %||% "", "")
M2$patch <- vapply(new[M2$id], function(r) paste(r$patch, collapse = "+"),
                   "")
cls <- vapply(seq_len(nrow(M2)), function(i) classify(M2[i, ]), "")
cat("classes changed by the new spell pass, against the lane's merged",
    "table:", sum(cls != M$class), "\n")
