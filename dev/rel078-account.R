# The counts of the merged build against the 0.77.0 release, file by
# file where 0.77.0's per-file counts can be had, and package by package,
# for one tier. Adapted from dev/rel077-account.R.
#
# 0.77.0's per-file logs went with its release tree; its package totals
# are in dev/round-20261028.md, from its tiers r077-plain (28966) and
# r077-ref (32103). Five test files changed after those tiers ran: the
# release review's fixes (core test-get-prior-route.R,
# test-prior-one-slot.R, test-arg-refusal.R; frmtmb.sample
# test-flat-prior.R) and ac4906bd (frmtmb.sample test-mo-coding.R).
# Their r077-ref counts are in that record (78, 82, 122, 121, 50). In
# the plain tier:
#   - the first four run with no skip without the gates (the before arm
#     below shows skip 0 and the gated count for each), so no block of
#     theirs is gated, and their r077-plain count is the r077-ref one;
#   - test-mo-coding.R has gated blocks (plain 9, skip 2), and ac4906bd
#     replaced one expectation by one (expect_identical() by
#     expect_equal(tolerance = 1e-12)), so its r077-plain count is
#     ac4906bd's in either tier.
# So:
#   - the test files the round changed and the five, as ac4906bd has
#     them, are rerun on rellib-r15 (r078-before077 gated, r078-before077p
#     plain; dev/rel078-log/before.jobs, from an export of ac4906bd in the
#     session scratchpad);
#   - the five carry 0.77.0's tier count to ac4906bd's ("0.77.0 after its
#     review" column);
#   - a file that is new has 0 before;
#   - every other file must account for nothing: each package's total
#     change must equal the sum of its listed files' deltas.
#   Rscript dev/rel078-account.R <tier: plain|ref> > dev/rel078-log/account-<tier>.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel078")
tier <- commandArgs(TRUE)[1]
stopifnot(tier %in% c("plain", "ref"))
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
run <- paste0("r078-", tier)
M <- rd(file.path("dev/ciharden-log", paste0(run, "-suite")))
B <- rd(file.path("dev/ciharden-log",
                  if (tier == "ref") "r078-before077-suite" else
                    "r078-before077p-suite"))
# 0.77.0's tier totals (dev/round-20261028.md)
P077 <- if (tier == "ref") {
  c(frmtmb = 22931, frmtmb.coupling = 557, frmtmb.eam = 1862,
    frmtmb.latent = 687, frmtmb.learn = 515, frmtmb.ode = 551,
    frmtmb.sample = 3887, frmtmb.spatial = 561, frmtmb.survival = 552)
} else {
  c(frmtmb = 20355, frmtmb.coupling = 551, frmtmb.eam = 1838,
    frmtmb.latent = 681, frmtmb.learn = 430, frmtmb.ode = 549,
    frmtmb.sample = 3532, frmtmb.spatial = 478, frmtmb.survival = 552)
}
stopifnot(sum(P077) == if (tier == "ref") 32103 else 28966)
F077 <- c(frmtmb = 285, frmtmb.coupling = 11, frmtmb.eam = 32,
          frmtmb.latent = 12, frmtmb.learn = 15, frmtmb.ode = 11,
          frmtmb.sample = 73, frmtmb.spatial = 8, frmtmb.survival = 15)
stopifnot(sum(F077) == 462)
key_b <- paste(B$pkg, B$file)
# r077's counts of the five files that changed after the tiers
R077 <- c("frmtmb test-get-prior-route" = 78,
          "frmtmb test-prior-one-slot" = 82,
          "frmtmb test-arg-refusal" = 122,
          "frmtmb.sample test-flat-prior" = 121,
          "frmtmb.sample test-mo-coding" = 50)
stopifnot(all(names(R077) %in% key_b))
post <- B[match(names(R077), key_b), ]
four <- names(R077) != "frmtmb.sample test-mo-coding"
stopifnot(all(post$skip[four] == 0L), all(post$fail == 0L))
if (tier == "plain") {
  R077[["frmtmb.sample test-mo-coding"]] <- post$pass[!four]
}
jobs <- read.table("dev/rel078-log/before.jobs", col.names = c("pkg", "path"))
cat("merged:", run, "files", nrow(M), "pass", sum(M$pass), "fail",
    sum(M$fail), "skip", sum(M$skip), "| 0.77.0: files", sum(F077),
    "pass", sum(P077), "\n\n")
cat("| file | 0.77.0 tier | 0.77.0 after its review | 0.78.0 |",
    "delta to 0.77.0 tier |\n|---|---|---|---|---|\n")
dl <- list()
exp_dir <- file.path("C:/Users/adf44/AppData/Local/Temp/12/claude",
  "c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105",
  "scratchpad/ac4906bd")
for (i in seq_len(nrow(M))) {
  key <- paste(M$pkg[i], M$file[i])
  j <- match(key, key_b)
  had <- file.exists(file.path(exp_dir,
    if (M$pkg[i] == "frmtmb") "tests/testthat" else
      file.path("extensions", M$pkg[i], "tests/testthat"),
    paste0(M$file[i], ".R")))
  if (had && is.na(j)) next  # unchanged since 0.77.0's tiers
  b077 <- if (!had) 0L else if (key %in% names(R077)) R077[[key]] else
    B$pass[j]
  bac <- if (!had) 0L else B$pass[j]
  cat("|", M$pkg[i], paste0("`", M$file[i], ".R`"), "|",
      if (!had) "new" else b077, "|", if (!had) "new" else bac, "|",
      M$pass[i], "|", M$pass[i] - b077, "|\n")
  dl[[M$pkg[i]]] <- c(dl[[M$pkg[i]]], M$pass[i] - b077)
}
cat("\n| package | files 0.77.0 | files 0.78.0 | pass 0.77.0 |",
    " pass 0.78.0 | change | listed files' sum | unexplained |\n",
    "|---|---|---|---|---|---|---|---|\n", sep = "")
for (p in names(P077)) {
  x <- M[M$pkg == p, ]
  s <- sum(dl[[p]])
  cat("|", p, "|", F077[[p]], "|", nrow(x), "|", P077[[p]], "|",
      sum(x$pass), "|", sum(x$pass) - P077[[p]], "|", s, "|",
      sum(x$pass) - P077[[p]] - s, "|\n")
}
cat("\n0.77.0 after its review, this tier:", sum(P077) +
      sum(post$pass - R077), "(the five files:",
    paste(sprintf("%s %+d", names(R077), post$pass - R077), collapse = ", "),
    ")\n")
cat("before arm: files", nrow(B), "fail", sum(B$fail), "skip", sum(B$skip),
    "| merged: fail", sum(M$fail), "skip", sum(M$skip), "\n")
