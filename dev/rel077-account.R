# Every-gate counts of the merged build against the 0.76.0 release, file
# by file where 0.76.0's per-file counts can be had, and package by
# package. 0.76.0's per-file logs went with its release tree; its
# package totals are in dev/round-20261027.md ("After the release
# review": r076-ref2, 31532 assertions). So:
#   - the test files the round changed, as 435e7ce6 has them, are rerun
#     on rellib-r14 (r077-before076, dev/rel077-log/before.jobs, from an
#     export of 435e7ce6's test directories);
#   - a file that is new has 0 before;
#   - every other file must account for nothing: each package's total
#     change must equal the sum of its changed files' deltas.
# The merged build is r077-ref.
#   Rscript dev/rel077-account.R > dev/rel077-log/account.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel077")
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
run <- Sys.getenv("REL077_ACCOUNT_RUN", "r077-ref")
M <- rd(file.path("dev/ciharden-log", paste0(run, "-suite")))
B <- rd("dev/ciharden-log/r077-before076-suite")
# test-arg-refusal.R skips two blocks without the package sources beside
# its tests; it was rerun from an export that has R/ and dev/ too
B2 <- rd("dev/ciharden-log/r077-before076b-suite")
B[match(paste(B2$pkg, B2$file), paste(B$pkg, B$file)),
  c("pass", "fail", "skip")] <- B2[, c("pass", "fail", "skip")]
# 0.76.0, every gate on, after its release fixes (dev/round-20261027.md)
P076 <- c(frmtmb = 22637, frmtmb.coupling = 557, frmtmb.eam = 1862,
          frmtmb.latent = 687, frmtmb.learn = 515, frmtmb.ode = 551,
          frmtmb.sample = 3610, frmtmb.spatial = 561, frmtmb.survival = 552)
stopifnot(sum(P076) == 31532)
F076 <- c(frmtmb = 281, frmtmb.coupling = 11, frmtmb.eam = 32,
          frmtmb.latent = 12, frmtmb.learn = 15, frmtmb.ode = 11,
          frmtmb.sample = 70, frmtmb.spatial = 8, frmtmb.survival = 15)
stopifnot(sum(F076) == 455)
jobs <- read.table("dev/rel077-log/before.jobs", col.names = c("pkg", "path"))
changed <- paste(jobs$pkg, sub("[.]R$", "", basename(jobs$path)))
cat("merged:", run, "files", nrow(M), "pass", sum(M$pass), "fail",
    sum(M$fail), "skip", sum(M$skip), "| 0.76.0: files", sum(F076), "pass",
    sum(P076), "\n\n")
cat("| file | 0.76.0 | 0.77.0 | delta |\n|---|---|---|---|\n")
dl <- list()
for (i in seq_len(nrow(M))) {
  key <- paste(M$pkg[i], M$file[i])
  j <- match(key, paste(B$pkg, B$file))
  isnew <- !key %in% changed && is.na(j) &&
    !file.exists(file.path("C:/Users/adf44/AppData/Local/Temp/12/claude",
      "c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105",
      "scratchpad/rel077-before076",
      if (M$pkg[i] == "frmtmb") "tests/testthat" else
        file.path("extensions", M$pkg[i], "tests/testthat"),
      paste0(M$file[i], ".R")))
  b <- if (isnew) 0L else if (!is.na(j)) B$pass[j] else NA
  if (is.na(b)) next
  cat("|", M$pkg[i], paste0("`", M$file[i], ".R`"), "|",
      if (isnew) "new" else b, "|", M$pass[i], "|", M$pass[i] - b, "|\n")
  dl[[M$pkg[i]]] <- c(dl[[M$pkg[i]]], M$pass[i] - b)
}
cat("\n| package | files 0.76.0 | files 0.77.0 | pass 0.76.0 | pass 0.77.0 |",
    "change | changed files' sum | unexplained |\n|---|---|---|---|---|---|---|---|\n")
for (p in names(P076)) {
  x <- M[M$pkg == p, ]
  s <- sum(dl[[p]])
  cat("|", p, "|", F076[[p]], "|", nrow(x), "|", P076[[p]], "|",
      sum(x$pass), "|", sum(x$pass) - P076[[p]], "|", s, "|",
      sum(x$pass) - P076[[p]] - s, "|\n")
}
cat("\nbefore arm: files", nrow(B), "fail", sum(B$fail), "| merged: fail",
    sum(M$fail), "skip", sum(M$skip), "\n")
