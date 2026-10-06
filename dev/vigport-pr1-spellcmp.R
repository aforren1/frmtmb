# Punch round 1: which spell-pass rows and counts moved when the two
# stale hypothesis patches (brms_overview.6.1, brms_distreg.5.1) left
# dev/brms-port/patches.R.
#
#   Rscript dev/vigport-pr1-spellcmp.R
#
# Compares dev/vigport-port-out/r5/results-merged-pre-pr1.rds (the spell
# pass with the patches) with results-merged.rds (without them).
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
o <- file.path(here, "vigport-port-out", "r5")
A <- readRDS(file.path(o, "results-merged-pre-pr1.rds"))
B <- readRDS(file.path(o, "results-merged.rds"))
stopifnot(identical(sort(A$id), sort(B$id)))
B <- B[match(A$id, B$id), ]
st <- A$status_spell != B$status_spell
cl <- A$class != B$class
cat("rows:", nrow(A), " spell status changed:", sum(st),
    " class changed:", sum(cl), "\n")
for (i in which(st | cl)) {
  cat(sprintf("%-26s %-5s spell %s -> %s; class %s -> %s\n", A$id[i],
              A$kind[i], A$status_spell[i], B$status_spell[i], A$class[i],
              B$class[i]))
}
for (k in c("model", "post", "other")) {
  # the headline convention: a SETUP-SKIP row is never run, so it is
  # not counted
  a <- A[A$kind == k & A$status_raw != "SETUP-SKIP", ]
  b <- B[B$kind == k & B$status_raw != "SETUP-SKIP", ]
  cat(sprintf("%-5s spell OK: %d -> %d of %d\n", k,
              sum(a$status_spell == "OK"), sum(b$status_spell == "OK"),
              nrow(a)))
}
