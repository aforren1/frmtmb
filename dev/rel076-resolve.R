# Consolidation edits of the 0.76.0 release tree, applied after the
# textual merge (dev/rel076-merge.sh) and before the NEWS merge
# (dev/rel076-news-merge.R). Each edit checks that its old text is
# there exactly once (or the count it names), so a rerun on a fresh
# merge either applies all of them or stops on the first that no longer
# fits. A file no lane changed is restored from its HEAD blob (70476dc3)
# before it is edited, so every edit here is applied once however often
# this runs. Every hand resolution of the round is here, each with its
# reason. Adapted from dev/rel075-resolve.R.
#
#   Rscript dev/rel076-resolve.R
rel <- Sys.getenv("REL076_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel076")
git_dir <- Sys.getenv("REL076_GIT", rel)
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
merged <- readLines("dev/rel076-log/merge-work/files")
copied <- unlist(lapply(list.files("dev/rel076-log", "^copied-.*[.]lst$",
                                   full.names = TRUE), readLines))
restore <- function(f) {
  if (f %in% c(merged, copied)) {
    stop(f, " is a merged or copied file; edit it, do not restore it")
  }
  from_head(f)
}
# Replace the one conflict block of `f` whose two sides contain `ours`
# and `theirs` (fixed strings) with `new` (a character vector of lines,
# or a function of the two sides' lines). Stops unless exactly one
# block matches.
resolve_block <- function(f, ours, theirs, new) {
  x <- rd(f)
  a <- which(startsWith(x, "<<<<<<< "))
  m <- which(x == "=======")
  z <- which(startsWith(x, ">>>>>>> "))
  stopifnot(length(a) == length(z))
  hit <- integer(0)
  for (i in seq_along(a)) {
    mid <- m[m > a[i] & m < z[i]]
    stopifnot(length(mid) == 1L)
    o <- paste(x[(a[i] + 1L):(mid - 1L)], collapse = "\n")
    t <- if (z[i] > mid + 1L) paste(x[(mid + 1L):(z[i] - 1L)],
                                    collapse = "\n") else ""
    if (grepl(ours, o, fixed = TRUE) && grepl(theirs, t, fixed = TRUE)) {
      hit <- c(hit, i)
    }
  }
  if (length(hit) != 1L) stop(f, ": expected 1 conflict block, found ",
                              length(hit), ": ", ours)
  i <- hit
  mid <- m[m > a[i] & m < z[i]]
  if (is.function(new)) {
    new <- new(x[seq_len(mid - a[i] - 1L) + a[i]],
               if (z[i] > mid + 1L) x[(mid + 1L):(z[i] - 1L)] else
                 character(0))
  }
  wr(c(x[seq_len(a[i] - 1L)], new, x[-seq_len(z[i])]), f)
  cat("resolved a conflict in", f, "\n")
}

# ---- 1. R/sampling-api.R: acprior x loglik --------------------------------
# Both lanes added names to the sampling-API page after
# rescor_row_loglik, acprior three (autocor_flat_priors(),
# autocor_bound_chart(), autocor_innov_sd()) and loglik one
# (autocor_row_loglik()), each in the @aliases list and in the list of
# the page's NULL block. The insertions share one anchor, so the merge
# could not place them; both are kept, acprior's first (the merge
# order). NAMESPACE and man/frmtmb-sampling-api.Rd are roxygen's.
both <- function(ours, theirs) c(ours, theirs)
resolve_block("R/sampling-api.R", "#' @aliases autocor_flat_priors",
              "#' @aliases autocor_row_loglik", both)
resolve_block("R/sampling-api.R", "#'   autocor_flat_priors,",
              "#'   autocor_row_loglik,", both)

# ---- 2. dev/upstream-bugs.md ----------------------------------------------
# mocode added brms-40 and brms-41; no other lane changed the file, so
# no number is taken twice. Checked, not edited.
ub <- rd("dev/upstream-bugs.md")
nums <- sub("^### (brms-[0-9]+)[.].*", "\\1",
            grep("^### brms-[0-9]+[.]", ub, value = TRUE))
if (anyDuplicated(nums)) stop("dev/upstream-bugs.md: a number twice: ",
                              nums[duplicated(nums)][1])
cat("dev/upstream-bugs.md:", length(nums), "entries, last",
    nums[length(nums)], ", no number twice\n")

# ---- 3. versions and floors -----------------------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
dcore <- "DESCRIPTION"
restore(dcore)
ver(dcore, "0.75.0", "0.76.0")
cm <- "codemeta.json"
restore(cm)
# no lane changed a dependency of core, so only the version moves
edit_file(cm, '"version": "0.75.0",', '"version": "0.76.0",')
# frmtmb.sample: calls autocor_flat_priors(), autocor_bound_chart(),
# autocor_innov_sd() and autocor_row_loglik(); minor (breaking default
# priors, log_lik() of residual covariance terms)
dsa <- "extensions/frmtmb.sample/DESCRIPTION"
restore(dsa)
ver(dsa, "0.22.0", "0.23.0")
edit_file(dsa, "frmtmb (>= 0.75.0)", "frmtmb (>= 0.76.0)")
# frmtmb.spatial: sar() registers with frmtmb_register_autocor(
# jacobian =), which 0.75.0 refuses as an unused argument at load
dsp <- "extensions/frmtmb.spatial/DESCRIPTION"
restore(dsp)
ver(dsp, "0.2.0", "0.3.0")
edit_file(dsp, "frmtmb (>= 0.75.0)", "frmtmb (>= 0.76.0)")

# ---- 4. the release harness reads rellib-r14 -------------------------------
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R", "dev/brmsport-bin2-all.sh",
            "dev/brmsport-bin2-record.sh", "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r13", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("rellib-r13", "rellib-r14", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}

# ---- 4b. integration tests the lanes could not write ----------------------
# cosy() + mo() through frm_sample() and log_lik(), and ar(cov = TRUE)
# under the default priors through log_lik(): the three changes meet
# only in the merged tree (dev/rel076-overlap-test.R holds the blocks,
# dev/rel076-overlap.R is the probe, which also covers sar()'s prior
# tables, a pair of packages no test file loads together)
blk <- "dev/rel076-overlap-test.R"
f <- "extensions/frmtmb.sample/tests/testthat/test-loglik-autocor.R"
stopifnot(f %in% copied)
x <- rd(f)
stopifnot(!any(grepl("the 0.76.0 consolidation", x, fixed = TRUE)))
wr(c(x, "", rd(blk)), f)
cat("appended", blk, "to", f, "\n")

# ---- 4c. the ported brms suite: brmsfit-methods:924 by policy -------------
# mocode's coding from the minimum makes the block "update has
# reasonable outputs" refuse its refit: its newdata hold Exp = 4 on
# every row, so mo(Exp) has one category. brms 2.23.0 builds that model
# only under testmode, which never samples: Stan rejects the data ("Jmo[1]
# is 0, but must be greater than or equal to 1"; dev/mocode-brms-const.R,
# dev/mocode-findings.md "A ported row that flips"). The review asked
# for a divergence by policy; the ledger (dev/brmsport-ledger.R) reads
# this file and stops on an assertion that fails without a verdict.
vm <- "dev/brmsport-verdicts-manual.tsv"
restore(vm)
x <- rd(vm)
stopifnot(!any(startsWith(x, "brmsfit-methods:924\t")))
wr(c(x, paste("brmsfit-methods:924", "frmtmb", "divergence", "policy",
              paste0("update(newdata =) refits with Exp = 4 on every row; ",
                     "mo() codes from the minimum as brms 2.23.0 does, so ",
                     "mo(Exp) has one category and frm() refuses it. brms ",
                     "builds the model only under testmode, which never ",
                     "samples: Stan rejects Jmo = 0 ",
                     "(dev/mocode-brms-const.R)"), sep = "\t")), vm)
cat("appended the verdict of brmsfit-methods:924 to", vm, "\n")

# ---- 5. release fixes -----------------------------------------------------
# the lanes' reviews and the coordinator's instructions; the fixes and
# their reasons are dev/rel076-release-fix.R, sourced here so a merge
# rerun applies them
source("dev/rel076-release-fix.R", local = TRUE)

# ---- 6. dev/test-backlog.md and dev/round-handoff.md ----------------------
bl <- "dev/test-backlog.md"
restore(bl)
bsec <- "dev/rel076-backlog-section.md"
if (file.exists(bsec)) {
  edit_file(bl, "\n## Reference\n",
            paste0("\n", paste(rd(bsec), collapse = "\n"),
                   "\n\n## Reference\n"))
} else cat("NOTE: no", bsec, "yet\n")
hf <- "dev/round-handoff.md"
restore(hf)
hsec <- "dev/rel076-handoff-section.md"
if (file.exists(hsec)) {
  edit_file(hf,
"0.74.0 release and 2026-10-24 at the 0.75.0 release. Read this, then",
"0.74.0 release, 2026-10-24 at the 0.75.0 release and 2026-10-27 at the
0.76.0 release. Read this, then")
  edit_file(hf, "\n## The 0.75.0 round, ",
            paste0("\n", paste(rd(hsec), collapse = "\n"),
                   "\n\n## The 0.75.0 round, "))
} else cat("NOTE: no", hsec, "yet\n")
cat("RESOLVE DONE\n")
