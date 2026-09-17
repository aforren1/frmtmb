## Punch round 2, MINOR 2: the fuzz tier (FRMTMB_FUZZ=true) generates
## binomial() rows without trials(). This runs exactly those rows of the
## tier's own plan (seed 20260901, size 300) through fuzz_run(), twice:
## with the helper as it now is (trials(nt) written, nt = 1), and with the
## formula builder reverted to the pre-fix text, so the fix is seen to
## matter on the lane build rather than assumed to.
## Usage (from the worktree root): Rscript dev/famlink-p2-fuzzbinom.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
suppressMessages(library(testthat))
env <- new.env(parent = asNamespace("frmtmb"))
sys.source("tests/testthat/helper-fuzz.R", envir = env)
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
rows <- plan[plan$family == "binomial" & plan$aterm != "trials", ]
cat("plan rows:", nrow(plan), "; binomial without a trials aterm:",
    nrow(rows), "; aterms:", paste(names(table(rows$aterm)),
                                   table(rows$aterm), collapse = ", "), "\n")
show <- function(label) {
  res <- env$fuzz_run(rows)
  sm <- env$fuzz_summary(res)
  cls <- vapply(sm$triaged, function(f) f$class, "")
  cat(sprintf("%-8s findings %d; by class: %s\n", label, length(cls),
              paste(names(table(cls)), table(cls), collapse = ", ")))
  new <- Filter(function(f) f$class %in% c("real_new", "generator"),
                sm$triaged)
  for (f in utils::head(new, 3)) {
    cat("   ", f$invariant, ":", substr(f$detail, 1, 100), "\n")
  }
}
cat("example call:", gsub("\n", " ", env$fuzz_call_text(as.list(rows[1, ]))),
    "\n")
show("fixed")
# the pre-fix builder: the aterm text without trials(nt)
env$fuzz_aterm_text <- function(sp) {
  fm <- env$fuzz_families[[sp$family]]
  if (identical(sp$aterm, "trunc")) {
    paste0("trunc(lb = ", format(fm$trunc_lb), ")")
  } else {
    env$fuzz_aterm[[sp$aterm]]$term
  }
}
cat("example call:", gsub("\n", " ", env$fuzz_call_text(as.list(rows[1, ]))),
    "\n")
show("prefix")
