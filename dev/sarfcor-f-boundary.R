# Fix round 1: dev/core-boundary.md's seam table. Part A's edit goes to
# the worktree and the stage-A copy, part B's to the worktree only.
# Rscript <this> <worktree> <stageA>
args <- commandArgs(TRUE)
edit <- function(root, pairs) {
  p <- file.path(root, "dev/core-boundary.md")
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
a <- list(c(
  "| `frmtmb_register_autocor()` (after 0.74.0) | a residual correlation term: one affine map from the response to independent innovations, such as `sar()` |",
  "| `frmtmb_register_autocor()` (after 0.74.0) | a residual correlation term: one affine map from the response to independent innovations, such as `sar()`; its `bounds` slot gives a bounded parameter brms's class prior and frmtmb.sample's flat default |
| `frmtmb_autocor_pars()` (after 0.74.0) | lists a fit's registered residual-term parameters with their intervals, for a package that would otherwise keep its own list of names (frmtmb.sample's `lagsar`, `errorsar`) |"))
b <- list(c(
  "| `frmtmb_register_covstruct(full_rank =)` (after 0.74.0) | whether a structure's covariance is positive definite at some value, which the weights check reads |",
  "| `frmtmb_register_covstruct(full_rank =)` (after 0.74.0) | whether a row's own effect can reach any value with a finite density, which the weights check reads; a structure that states nothing leaves the check open |"))
for (root in args) edit(root, a)
edit(args[1], b)
cat("ok\n")
