p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/R/weights-unbounded.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("# Covariance structures with a positive-definite covariance at some
# parameter value, so the effects can take any value with a finite
# density. rr() and the structures with a fixed or possibly singular
# covariance are left out.
unbounded_full_rank_blocks <- c(
  \"us\", \"diag\", \"homdiag\", \"cs\", \"homcs\", \"ar1\", \"hetar1\", \"ou\", \"toep\",
  \"homtoep\", \"exp\", \"gau\", \"mat\")
",
"# Which covariance structures can give a row its own effect at any value
# with a finite density is a property of the structure, read through
# cs_full_rank() (R/covstruct.R): a built-in entry states it, and a
# registered structure states it or has it read off its covariance.
")
rep1("  owned_blocks <- Filter(function(bk) {
    bk[[\"covstruct\"]] %in% unbounded_full_rank_blocks
  }, frame[[\"re_blocks\"]])",
"  owned_blocks <- Filter(cs_full_rank, frame[[\"re_blocks\"]])")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
