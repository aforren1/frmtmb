# Consolidation edits of the 0.78.0 release tree, applied after the
# textual merge (dev/rel078-merge.sh) and before the bin-2 generator and
# the NEWS merge (dev/rel078-news-merge.R). Each edit checks that its old
# text is there exactly once (or the count it names), so a rerun on a
# fresh merge either applies all of them or stops on the first that no
# longer fits. A file no lane changed is restored from its HEAD blob
# (ac4906bd) before it is edited, so every edit here is applied once
# however often this runs. Every hand resolution of the round is here,
# each with its reason. Adapted from dev/rel077-resolve.R.
#
#   Rscript dev/rel078-resolve.R
rel <- Sys.getenv("REL078_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel078")
git_dir <- Sys.getenv("REL078_GIT", rel)
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
setwd(rel)
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)),
               stdout = TRUE)
  wr(sub("\r$", "", x), f)
}
edit_file <- function(f, old, new, n = 1L) {
  s <- paste(rd(f), collapse = "\n")
  hits <- gregexpr(old, s, fixed = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ",
                   substr(old, 1, 70))
  s <- gsub(old, new, s, fixed = TRUE)
  wr(strsplit(s, "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}
merged <- readLines("dev/rel078-log/merge-work/files")
copied <- unlist(lapply(list.files("dev/rel078-log", "^copied-.*[.]lst$",
                                   full.names = TRUE), readLines))
restore <- function(f) {
  if (f %in% c(merged, copied)) {
    stop(f, " is a merged or copied file; edit it, do not restore it")
  }
  from_head(f)
}

# ---- 1. the textual conflicts -------------------------------------------
# Only the two NEWS files conflict: meme and rescorpp each insert a
# "(development version)" section at line 1 of core NEWS.md and of
# frmtmb.sample's NEWS.md. dev/rel078-news-merge.R rebuilds both top
# sections from the lanes' own sections (one heading, one "Needs" lead),
# and takes the released sections below from the HEAD blob, so the
# conflict blocks never reach the result. No code file is changed by two
# lanes (dev/rel078-log/overlap-files.txt): checked here.
st <- lapply(c(equiv = "equiv", meme = "meme", rescorpp = "rescorpp"),
             function(l) {
  f <- file.path("dev/rel078-log/merge-work", paste0(l, ".status"))
  if (!file.exists(f) || !file.size(f)) return(character(0))
  sub("^[A-Z]\t", "", readLines(f))
})
both <- intersect(st$meme, st$rescorpp)
stopifnot(setequal(both,
                    c("NEWS.md", "extensions/frmtmb.sample/NEWS.md")),
          length(st$equiv) == 0L)
wr(c(paste("changed by meme:", length(st$meme)),
     paste("changed by rescorpp:", length(st$rescorpp)),
     paste("changed by equiv:", length(st$equiv)),
     paste("changed by meme and rescorpp:", paste(both, collapse = ", "))),
   "dev/rel078-log/overlap-files.txt")
cat("tracked files two lanes change:", paste(both, collapse = ", "), "\n")

# ---- 2. dev/upstream-bugs.md: the two brms defects rescorpp found ------
# not posted; no lane changed the file; the numbers stay unique
ub <- "dev/upstream-bugs.md"
restore(ub)
edit_file(ub, "\n## RTMB\n", paste0("\n", paste(rd(
  "dev/rel078-upstream-brms43-44.md"), collapse = "\n"), "\n\n## RTMB\n"))
nums <- sub("^### (brms-[0-9]+)[.].*", "\\1",
            grep("^### brms-[0-9]+[.]", rd(ub), value = TRUE))
if (anyDuplicated(nums)) stop("dev/upstream-bugs.md: a number twice: ",
                              nums[duplicated(nums)][1])
stopifnot(identical(tail(nums, 2L), c("brms-43", "brms-44")))
cat("dev/upstream-bugs.md:", length(nums), "entries, last",
    nums[length(nums)], ", no number twice\n")

# ---- 3. versions and floors ---------------------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
dcore <- "DESCRIPTION"
restore(dcore)
# minor, breaking: brms's me() prior classes, the sampling-route group
# rows, a new sampling-API export
ver(dcore, "0.77.0", "0.78.0")
cm <- "codemeta.json"
restore(cm)
# no lane changed a dependency of core, so only the version moves
edit_file(cm, '"version": "0.77.0",', '"version": "0.78.0",')
# frmtmb.sample: calls predict_rescor_draw() and rescor_row_loglik(resp
# =), which 0.77.0 lacks, and writes class "sdme" and "corme" priors,
# which 0.77.0's set_prior() refuses; minor, breaking. Measured by
# dev/rel078-floor.sh
dsa <- "extensions/frmtmb.sample/DESCRIPTION"
restore(dsa)
ver(dsa, "0.24.0", "0.25.0")
edit_file(dsa, "frmtmb (>= 0.77.0)", "frmtmb (>= 0.78.0)")
# the other seven: dev/rel078-versions.R, which the measurement of their
# tests on the old and the new core decides (dev/round-20261029.md)
vf <- "dev/rel078-versions.R"
if (file.exists(vf)) source(vf, local = TRUE)

# ---- 4. the release harness reads rellib-r16 -----------------------------
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R", "dev/brmsport-bin2-all.sh",
            "dev/brmsport-bin2-record.sh", "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r15", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("rellib-r15", "rellib-r16", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}
# Rtools last on PATH in the per-file and the one-session runners: with
# its directories first, an R child segfaulted at exit in the 0.77.0
# bin-2 pipeline (dev/round-handoff.md, "What the 0.77.0 round is
# evidence for"); Git's bash and coreutils then stay first, and R still
# finds Rtools' make and compilers further down PATH
rss <- "dev/release/run-sessions.sh"
if (!rss %in% merged) restore(rss)
for (f in c("dev/ciharden-scan.sh", rss)) {
  edit_file(f,
            "export PATH=\"$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH\"",
            "export PATH=\"$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin\"")
}

# ---- 5. release fixes -----------------------------------------------------
# meme's review (S3) and the coordinator's instructions; the fixes and
# their reasons are dev/rel078-release-fix.R, sourced here so a merge
# rerun applies them
source("dev/rel078-release-fix.R", local = TRUE)

# ---- 6. dev/test-backlog.md and dev/round-handoff.md ---------------------
bl <- "dev/test-backlog.md"
restore(bl)
bsec <- "dev/rel078-backlog-section.md"
if (file.exists(bsec)) {
  edit_file(bl, "\n## Reference\n",
            paste0("\n", paste(rd(bsec), collapse = "\n"),
                   "\n\n## Reference\n"))
} else cat("NOTE: no", bsec, "yet\n")
hf <- "dev/round-handoff.md"
restore(hf)
hsec <- "dev/rel078-handoff-section.md"
if (file.exists(hsec)) {
  edit_file(hf,
"0.76.0 release and 2026-10-28 at the 0.77.0 release. Read this, then",
"0.76.0 release, 2026-10-28 at the 0.77.0 release and 2026-10-29 at the
0.78.0 release. Read this, then")
  edit_file(hf, "\n## The 0.77.0 round, ",
            paste0("\n", paste(rd(hsec), collapse = "\n"),
                   "\n\n## The 0.77.0 round, "))
} else cat("NOTE: no", hsec, "yet\n")
cat("RESOLVE DONE\n")
