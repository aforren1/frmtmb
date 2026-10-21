# Paste the generated blocks of dev/cdfs-findings.md from the logs, so
# that no number in it is typed. Each block sits between
#   <!-- GEN:<name> -->  and  <!-- /GEN:<name> -->
# and is replaced in place.
#
#   Rscript dev/cdfs-findings-build.R
f <- "dev/cdfs-findings.md"
md <- readLines(f)
put <- function(md, name, lines) {
  a <- grep(sprintf("^<!-- GEN:%s -->$", name), md)
  b <- grep(sprintf("^<!-- /GEN:%s -->$", name), md)
  if (length(a) != 1L || length(b) != 1L) {
    warning("no block ", name)
    return(md)
  }
  c(md[seq_len(a)], lines, md[b:length(md)])
}
rd <- function(p) if (file.exists(p)) readLines(p, warn = FALSE) else
  paste("MISSING:", p)
fence <- function(x) c("```", x, "```")

# brms's survey: one row per family and term
sv <- rd("dev/cdfs-out/survey.txt")
hd <- grep("^== ", sv)
rows <- character()
for (i in seq_along(hd)) {
  h <- trimws(sub(" [|] ", ", ", sub("^== ", "", sv[hd[i]])))
  body <- sv[(hd[i] + 1):(if (i < length(hd)) hd[i + 1] - 1 else length(sv))]
  err <- grep("ERROR", body, value = TRUE)
  fun <- regmatches(body, regexpr("[a-z_0-9]+_l(c)?cdf\\([^|]*", body))
  fun <- unique(sub("[(].*", "", fun))
  st <- if (length(err)) trimws(sub(".*ERROR: ", "refused: ", err[1])) else
    paste(fun, collapse = ", ")
  rows <- c(rows, sprintf("| %s | %s |", h, st))
}
md <- put(md, "survey", c("| family and term | brms 2.23.0 |", "|---|---|", rows))

md <- put(md, "probe", fence(c(
  "-- dev/cdfs-out/probe1.txt (which RTMB CDFs tape; gradient vs FD)",
  rd("dev/cdfs-out/probe1.txt"),
  "-- dev/cdfs-out/probe3.txt (tails; first lines)",
  utils::head(rd("dev/cdfs-out/probe3.txt"), 40))))

md <- put(md, "acc", fence(rd("dev/cdfs-out/acc.txt")))
md <- put(md, "grad", fence(rd("dev/cdfs-out/grad.txt")))
md <- put(md, "brmslik", fence(grep("^lib|^RESULT|^DETAIL|^LPCHECK|^CASE",
                                    rd("dev/cdfs-out/brms-lik-report.txt"),
                                    value = TRUE)))
md <- put(md, "flips", fence(rd("dev/cdfs-out/bin2-flips.txt")))
md <- put(md, "check", fence(grep("^Status|^[*] (checking|DONE).*(NOTE|WARNING|ERROR)",
  rd("dev/cdfs-check/check-frmtmb.txt"), value = TRUE)))
md <- put(md, "fuzzbase", fence(c(
  "-- dev/cdfs-out/fuzz-base-r10.txt (rellib-r10)",
  grep("^lib|^specs|^real_new|^    spec|errored|failed",
       rd("dev/cdfs-out/fuzz-base-r10.txt"), value = TRUE),
  "-- dev/cdfs-out/fuzz-base-mine.txt (wt-cdfs-lib)",
  grep("^lib|^specs|^real_new|^    spec|errored|failed",
       rd("dev/cdfs-out/fuzz-base-mine.txt"), value = TRUE))))

# test results, one line per file, from the per-file logs
res <- function(dir) {
  fs <- list.files(dir, "[.]txt$", full.names = TRUE)
  out <- vapply(fs, function(x) {
    r <- grep("^RESULT", readLines(x, warn = FALSE), value = TRUE)
    if (length(r)) utils::tail(r, 1) else paste("NO RESULT", basename(x))
  }, "")
  unname(out)
}
tot <- function(lines) {
  num <- function(k) sum(as.numeric(sub(".*=", "", regmatches(lines,
    regexpr(paste0(k, "=[0-9]+"), lines)))))
  sprintf("files %d: pass %d fail %d err %d skip %d warn %d",
          length(lines), num("pass"), num("fail"), num("err"), num("skip"),
          num("warn"))
}
# the fix round's runs; the first build's are kept in the logs named in
# the first report and are not repeated here
tiers <- c(ungated = "dev/release/cdfs-fx-all-files",
           # every file of the eight suites with every gate on, under the
           # reference BLAS: the gated tier and the reference arm at once
           gated_reference = "dev/ciharden-log/cdfs-fx-ref-suite",
           openblas326 = "dev/ciharden-log/cdfs-fx-ob326-suite",
           openblas332 = "dev/ciharden-log/cdfs-fx-ob332-suite",
           # test-cens-trunc.R after its five-row trunc() fit, which warned
           # under OpenBLAS, was given forty rows
           cens_trunc_ob326 = "dev/ciharden-log/cdfs-fx-ct-ob0326-suite",
           cens_trunc_ob332 = "dev/ciharden-log/cdfs-fx-ct-ob0332-suite",
           cens_trunc_ref = "dev/ciharden-log/cdfs-fx-ct-ref-suite",
           # the new gated test on rellib-r10: it must fail there first
           base_r10_brmslik = "dev/release/cdfs-base-brmslik-files")
blk <- character()
for (nm in names(tiers)) {
  if (!dir.exists(tiers[[nm]])) {
    blk <- c(blk, paste(nm, ": not run"))
    next
  }
  r <- res(tiers[[nm]])
  bad <- r[!grepl("fail=0 err=0", r)]
  blk <- c(blk, paste0(nm, " (", tiers[[nm]], "): ", tot(r)), bad)
}
blk <- c(blk, "test-cdfs.R on rellib-r10 (dev/cdfs-out/base-test-cdfs.txt):",
         grep("^RESULT", rd("dev/cdfs-out/base-test-cdfs.txt"), value = TRUE))
md <- put(md, "tests", fence(blk))
md <- put(md, "sessions", fence(grep("SESSION|LEAK|ran ",
  rd("dev/release/cdfs-fx-sess.log"), value = TRUE)))

# fix round 1
md <- put(md, "fix1fuzz", fence(c(
  grep("^RESULT", rd("dev/ciharden-log/cdfs-fx-ref-suite/frmtmb--test-fuzz.txt"),
       value = TRUE),
  "-- dev/cdfs-fuzz-base.R on this build, the three specs of section 8",
  grep("^specs|^real_new|^known|^    spec|errored|failed",
       rd("dev/cdfs-out/fuzz-base-mine.txt"), value = TRUE))))
md <- put(md, "fix1cost", fence(c(rd("dev/cdfs-out/fix1-cost.txt"),
  "-- rellib-r10", rd("dev/cdfs-out/fix1-imp-rellib-r10.txt"),
  "-- this build", rd("dev/cdfs-out/fix1-imp-wt-cdfs-lib.txt"))))
md <- put(md, "fix1igamma", fence(rd("dev/cdfs-out/fix1-igamma.txt")))
md <- put(md, "fix1conv", fence(rd("dev/cdfs-out/fix1-conv.txt")))
md <- put(md, "fix1exg", fence(c("-- rellib-r10",
  grep("^u |^seed|^lib", rd("dev/cdfs-out/fix1-exgfit-r10.txt"), value = TRUE),
  "-- this build",
  grep("^u |^seed|^lib", rd("dev/cdfs-out/fix1-exgfit-lane.txt"), value = TRUE))))
acc <- rd("dev/cdfs-out/rerev/acc.txt")
w <- suppressWarnings(as.numeric(sub(".* worst ([^ ]+) .*", "\\1",
  grep(" worst ", acc, value = TRUE))))
md <- put(md, "fix1rerev", fence(c(
  "-- acc.txt (the review's tail grid): comparisons, and at most 1e-12",
  sprintf("%d comparisons, %d at most 1e-12; the rest:", length(w),
          sum(w <= 1e-12)),
  grep(" worst ", acc, value = TRUE)[w > 1e-12],
  "-- zitrunc4.txt (finding 2)",
  grep("diff|hand", rd("dev/cdfs-out/rerev/zitrunc4.txt"), value = TRUE),
  "-- trlb.txt (finding 4)", rd("dev/cdfs-out/rerev/trlb.txt"),
  "-- fuzzweib.txt (finding 4, the unpatched spec)",
  grep("as is|fits:", rd("dev/cdfs-out/rerev/fuzzweib.txt"), value = TRUE),
  "-- grad5.txt (finding 3)",
  grep("codes kept", rd("dev/cdfs-out/rerev/grad5.txt"), value = TRUE),
  "-- grad4.txt (finding 5)",
  grep("worst", rd("dev/cdfs-out/rerev/grad4.txt"), value = TRUE),
  "-- exgdens2.txt (finding 5)",
  grep("sigma / beta|directly", rd("dev/cdfs-out/rerev/exgdens2.txt"),
       value = TRUE))))
writeLines(md, f)
cat("rebuilt", f, "\n")
# which blocks of the new and changed tests fail on the pre-fix build
pre <- unlist(lapply(list.files("dev/release/cdfs-pre-files", "[.]txt$",
                                full.names = TRUE), readLines, warn = FALSE))
pre <- grep("^DETAIL (FAIL|ERROR)", pre, value = TRUE)
blocks <- unique(sub("^DETAIL (FAIL|ERROR) [[]([^]]*)[]].*", "\\2", pre))
md <- readLines(f)
md <- put(md, "fix1pre", fence(c(grep("RESULT", rd("dev/release/cdfs-pre.log"),
                                      value = TRUE),
  utils::head(rd("dev/release/cdfs-pre-files/frmtmb--test-cdfs.txt"), 1),
  sprintf("%d blocks fail or error there:", length(blocks)), sort(blocks))))
writeLines(md, f)
# fix round 2
pre2 <- unlist(lapply(list.files("dev/release/cdfs-pre2-files", "[.]txt$",
                                 full.names = TRUE), readLines, warn = FALSE))
pre2 <- grep("^DETAIL (FAIL|ERROR)", pre2, value = TRUE)
b2 <- unique(sub("^DETAIL (FAIL|ERROR) [[]([^]]*)[]].*", "\\2", pre2))
keep <- function(p, re) grep(re, rd(p), value = TRUE)
md <- readLines(f)
md <- put(md, "fix2", fence(c(
  "-- the round's new tests on the pre-fix build (wt-cdfs-rev-lib)",
  grep("RESULT", rd("dev/release/cdfs-pre2.log"), value = TRUE),
  sprintf("%d blocks fail or error there:", length(b2)), sort(b2),
  "-- N1: the review's sumtail.R and sumfit.R on this build",
  keep("dev/cdfs-out/rerev/sumtail.txt", "interval|got"),
  keep("dev/cdfs-out/rerev/sumfit.txt", "intercept|independent|frmtmb"),
  "-- N1: the review's compois.R on this build",
  keep("dev/cdfs-out/rerev/compois.txt", "^codes"),
  "-- N3: dev/cdfs-fix2-simcost.R",
  keep("dev/cdfs-out/fix2-simcost.txt", "inv_after"),
  "-- N4: dev/cdfs-fix1-imp.R, this build and rellib-r10",
  keep("dev/cdfs-out/fix1-imp-wt-cdfs-lib.txt", "importance|peak"),
  keep("dev/cdfs-out/fix1-imp-rellib-r10.txt", "importance|peak"),
  "-- N5: dev/cdfs-fix2-brms-top.R (brms 2.23.0)",
  keep("dev/cdfs-out/fix2-brms-top.txt", "brm ran|brm ERROR"))))
writeLines(md, f)
md <- readLines(f)
ft <- rd("dev/cdfs-out/fix2-fuzztrunc.txt")
md <- put(md, "fix2fuzz", fence(c(
  grep("RESULT", rd("dev/release/cdfs-g6.log"), value = TRUE),
  grep("REAL-NEW|spec:", unlist(strsplit(paste(grep("^DETAIL",
    unlist(lapply(list.files("dev/release/cdfs-g6-files", full.names = TRUE),
                  readLines, warn = FALSE)), value = TRUE), collapse = " | "),
    " [|] ")), value = TRUE),
  "-- dev/cdfs-fix2-fuzztrunc.R, before the entry: classes, and the",
  "   REAL-NEW records' messages",
  ft[grep("^cls", ft) + 0:2],
  unique(grep("errored", ft, value = TRUE)))))
writeLines(md, f)
