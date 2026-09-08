# An extension that depends on another extension needs its check job to
# INSTALL that sibling from the checkout, and to re-run when it changes.
#
# WHY THIS EXISTS. A sibling is in no repository. Every extension job
# resolves the rest of its dependencies with `install.packages()`, which
# cannot find one, so the job dies at
#
#   ERROR: dependency 'frmtmb.eam' is not available for package
#   'frmtmb.learn'
#
# before a single test runs. That is not hypothetical: frmtmb.learn 0.2.0
# took `frmtmb.eam` into Imports and its job failed exactly that way,
# and while a lane was adding two more siblings to frmtmb.sample's
# Suggests its reviewer measured the same ERROR again. Both were caught
# by a human reading a red run. This file catches it in the suite.
#
# The path filter half matters for a different reason: without it a
# change to the sibling does not trigger the dependent's job at all, so
# the pair goes untested rather than failing loudly.
#
# WHAT THIS DOES NOT CHECK. Whether the sibling's own dependencies are
# installed before it, and whether the install order is right. Both are
# real and both need a YAML parse rather than a grep, which is not worth
# a dependency here; the two jobs that install a sibling do it after
# core and that is asserted by nothing but review.

ci_root <- function() {
  root <- testthat::test_path("..", "..")
  ok <- file.exists(file.path(root, "DESCRIPTION")) &&
    any(trimws(readLines(file.path(root, "DESCRIPTION"), n = 5L)) ==
          "Package: frmtmb") &&
    dir.exists(file.path(root, "extensions")) &&
    dir.exists(file.path(root, ".github", "workflows"))
  if (ok) root else NULL
}

# The workflow file for an extension is named from the package with dots
# turned into dashes. frmtmb.eam's job is the exception: it kept the
# name it had before the package was renamed from frmtmb.ddm.
ci_workflow_of <- function(root, pkg) {
  cand <- file.path(root, ".github", "workflows",
                    paste0("check-", gsub(".", "-", pkg, fixed = TRUE),
                           ".yaml"))
  if (file.exists(cand)) return(cand)
  alt <- list.files(file.path(root, ".github", "workflows"),
                    pattern = "^check-frmtmb", full.names = TRUE)
  hit <- alt[vapply(alt, function(f) {
    any(grepl(paste0("name: ", pkg, "$"), readLines(f))) }, logical(1))]
  if (length(hit) == 1L) hit else NA_character_
}

ci_siblings_of <- function(dir, pkg) {
  d <- read.dcf(file.path(dir, "DESCRIPTION"))
  f <- intersect(c("Depends", "Imports", "Suggests"), colnames(d))
  if (!length(f)) return(character(0))
  p <- unlist(strsplit(paste(d[, f], collapse = ","), ","))
  p <- trimws(sub("[(].*", "", p))
  p <- p[grepl("^frmtmb[.]", p)]
  setdiff(unique(p[nzchar(p)]), pkg)
}

test_that("an extension's check job installs every sibling it depends on", {
  root <- ci_root()
  skip_if_not(!is.null(root), "not a monorepo checkout")

  dirs <- list.dirs(file.path(root, "extensions"), recursive = FALSE)
  dirs <- dirs[file.exists(file.path(dirs, "DESCRIPTION"))]
  skip_if_not(length(dirs) > 0L, "no extensions in this tree")

  checked <- 0L
  for (d in dirs) {
    pkg <- basename(d)
    sibs <- ci_siblings_of(d, pkg)
    if (!length(sibs)) next
    wf <- ci_workflow_of(root, pkg)
    expect_false(is.na(wf),
                 info = paste0(pkg, " depends on ",
                               paste(sibs, collapse = ", "),
                               " but has no check workflow"))
    if (is.na(wf)) next
    y <- readLines(wf)
    for (s in sibs) {
      checked <- checked + 1L
      expect_true(
        any(grepl(paste0("R CMD INSTALL.*extensions/", s, "\\b"), y)),
        info = paste0(pkg, " depends on ", s, ", which is in no ",
                      "repository, but ", basename(wf), " never installs ",
                      "it from the checkout. `install.packages()` cannot ",
                      "resolve it and the job will fail at `dependency '",
                      s, "' is not available` before any test runs. See ",
                      "check-frmtmb-learn.yaml for the pattern."))
      expect_true(
        any(grepl(paste0("'extensions/", s, "/[*][*]'"), y, fixed = FALSE)),
        info = paste0(basename(wf), " does not list 'extensions/", s,
                      "/**' in its path filters, so a change to ", s,
                      " will not re-run ", pkg, "'s check even though ",
                      pkg, " depends on it."))
    }
  }
  # the assertions above are all inside a conditional, so prove the loop
  # reached them rather than passing on an empty sweep
  expect_gt(checked, 0L)
})
