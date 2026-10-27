# Every-gate counts of the merged build against the 0.75.0 release, file
# by file where 0.75.0's per-file logs can be had, and package by
# package. 0.75.0's per-file logs went with its release tree; its
# package totals are in dev/round-20261024.md ("After the release
# review": r075-ref plus the three files the release fixes changed,
# 30912 assertions). So:
#   - the files the round changed, as 70476dc3 has them, are rerun on
#     rellib-r13 (r076-before075, dev/rel076-log/before.jobs);
#   - a file that is new has 0 before;
#   - every other file must account for nothing: each package's total
#     change must equal the sum of its changed files' deltas.
# The merged build is r076-ref, with the two brms-suite files the
# ledger regenerated read from their rerun (r076-regen).
#   Rscript dev/rel076-account.R > dev/rel076-log/account.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel076")
rd <- function(dir) {
  fs <- list.files(dir, "[.]txt$", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    r <- tail(grep("^RESULT ", x, value = TRUE), 1)
    g <- function(k) {
      if (!length(r)) return(NA_integer_)
      as.integer(sub(paste0(".*", k, "=([0-9]+).*"), "\\1", r))
    }
    b <- sub("[.]txt$", "", basename(f))
    data.frame(pkg = sub("--.*$", "", b), file = sub("^.*--", "", b),
               pass = g("pass"), fail = g("fail"), skip = g("skip"))
  }))
}
M <- rd("dev/ciharden-log/r076-ref-suite")
G <- rd("dev/ciharden-log/r076-regen-suite")
k <- paste(M$pkg, M$file)
M[match(paste(G$pkg, G$file), k), c("pass", "fail", "skip")] <-
  G[, c("pass", "fail", "skip")]
B <- rd("dev/ciharden-log/r076-before075-suite")
# 0.75.0, every gate on, after its release fixes (dev/round-20261024.md)
P075 <- c(frmtmb = 22331, frmtmb.coupling = 557, frmtmb.eam = 1862,
          frmtmb.latent = 687, frmtmb.learn = 515, frmtmb.ode = 551,
          frmtmb.sample = 3329, frmtmb.spatial = 528, frmtmb.survival = 552)
stopifnot(sum(P075) == 30912)
F075 <- c(frmtmb = 279, frmtmb.coupling = 11, frmtmb.eam = 32,
          frmtmb.latent = 12, frmtmb.learn = 15, frmtmb.ode = 11,
          frmtmb.sample = 67, frmtmb.spatial = 8, frmtmb.survival = 15)
new <- c("frmtmb test-mo-coding", "frmtmb test-autocor-loglik",
         "frmtmb.sample test-mo-coding",
         "frmtmb.sample test-default-priors-autocor",
         "frmtmb.sample test-loglik-autocor")
cat("merged: files", nrow(M), "pass", sum(M$pass), "fail", sum(M$fail),
    "skip", sum(M$skip), "| 0.75.0: files", sum(F075), "pass",
    sum(P075), "\n\n")
cat("| file | 0.75.0 | 0.76.0 | delta |\n|---|---|---|---|\n")
dl <- list()
for (i in seq_len(nrow(M))) {
  key <- paste(M$pkg[i], M$file[i])
  j <- match(key, paste(B$pkg, B$file))
  b <- if (key %in% new) 0L else if (!is.na(j)) B$pass[j] else NA
  if (is.na(b)) next
  cat("|", M$pkg[i], paste0("`", M$file[i], ".R`"), "|",
      if (key %in% new) "new" else b, "|", M$pass[i], "|",
      M$pass[i] - b, "|\n")
  dl[[M$pkg[i]]] <- c(dl[[M$pkg[i]]], M$pass[i] - b)
}
cat("\n| package | files 0.75.0 | files 0.76.0 | pass 0.75.0 | pass 0.76.0 |",
    "change | changed files' sum | unexplained |\n|---|---|---|---|---|---|---|---|\n")
for (p in names(P075)) {
  x <- M[M$pkg == p, ]
  s <- sum(dl[[p]])
  cat("|", p, "|", F075[[p]], "|", nrow(x), "|", P075[[p]], "|",
      sum(x$pass), "|", sum(x$pass) - P075[[p]], "|", s, "|",
      sum(x$pass) - P075[[p]] - s, "|\n")
}
cat("\nfail", sum(M$fail), "skip", sum(M$skip), "\n")
