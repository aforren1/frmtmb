# The core/extensions boundary, at zero.
#
# This test began as a growth-only ratchet: dev/structured-family-
# protocol.md replaced the per-family branches in core with one
# `fam$structure` slot, and the pinned inventory measured the distance
# left to go, shrinking release by release (25 pairs at v0.44.0, 15 at
# v0.45.0, 4 at v0.46.0). The extraction round emptied it: the ODE seam
# left for extensions/frmtmb.ode, the sampling surface for
# extensions/frmtmb.sample, and with them the last tokens. The test is
# now the boundary itself: NO core file may name hmm, lca, their frame
# slots, or the ODE symbols, outside the two family homes that remain
# in core until protocol step 10 moves them out.
#
# `mixture` is deliberately NOT policed. It stays in core as the
# protocol's reference implementation of the structured-family seam.

test_that("no core file carries a structured-family or ODE reference", {
  # positive identification of the package SOURCE tree, not just an R/
  # directory: one CI layout offered an existing-but-empty ../../R, and
  # the guard must fail closed (skip) rather than open (assert nothing)
  rdir <- testthat::test_path("..", "..", "R")
  desc <- testthat::test_path("..", "..", "DESCRIPTION")
  is_src <- file.exists(desc) &&
    any(trimws(readLines(desc, n = 5L)) == "Package: frmtmb") &&
    dir.exists(rdir) &&
    file.exists(file.path(rdir, "objective.R"))
  skip_if_not(is_src,
              "package sources are not available (installed-package run)")

  # each token maps to the files still allowed to spell it: the two
  # family homes that live in core until step 10. R/ode.R is gone, so
  # the ODE symbols have no allowed home at all.
  tokens <- list(
    hmm     = c("hmm.R", "lca.R"),
    lca     = c("hmm.R", "lca.R"),
    mix_g   = c("hmm.R", "lca.R"),
    hmm_g   = c("hmm.R", "lca.R"),
    frm_ode = character(0),
    RTMBode = character(0)
  )

  # R's own \b counts `_` as a word character, so `\bhmm\b` would miss
  # hmm_check_fit() - exactly the symbols this test exists to see. The
  # boundary is therefore spelled against letters and digits only, which
  # still keeps `lca` out of `allcaps` and `hmm` out of a longer word.
  boundary_scan <- function(rdir, exempt = TRUE) {
    files <- sort(basename(list.files(rdir, pattern = "\\.R$")))
    out <- integer(0)
    for (tok in names(tokens)) {
      pat <- paste0("(?<![A-Za-z0-9])", tok, "(?![A-Za-z0-9])")
      scan_files <- if (exempt) setdiff(files, tokens[[tok]]) else files
      for (f in scan_files) {
        ln <- readLines(file.path(rdir, f), warn = FALSE)
        # roxygen and comment lines are prose about the seam, not the
        # seam itself. A trailing comment on a code line still counts.
        ln <- ln[!grepl("^\\s*#", ln)]
        n <- sum(vapply(gregexpr(pat, ln, perl = TRUE),
                        function(m) if (m[1L] == -1L) 0L else length(m),
                        0L))
        if (n > 0L) out[paste0(f, "|", tok)] <- n
      }
    }
    if (length(out)) out[order(names(out))] else out
  }

  # positive control: a scanner that matched nothing anywhere would
  # pass the boundary vacuously, so prove it still sees the tokens in
  # the exempt family homes before trusting the empty result
  control <- boundary_scan(rdir, exempt = FALSE)
  expect_gt(sum(control), 0L)

  found <- boundary_scan(rdir)
  hits <- sprintf("  %s: %d", names(found), found)
  expect_identical(
    names(found), NULL,
    info = paste0(
      "A core file names a structured family or an ODE symbol:\n",
      paste(hits, collapse = "\n"),
      "\n\nThe boundary is zero: put the branch behind the family's ",
      "own fam$structure slot, register it through the exported seams ",
      "(frmtmb_register_frame_check, frmtmb_register_compat), or move ",
      "the code into R/hmm.R, R/lca.R or the extension packages. ",
      "mixture() is not policed - it stays in core as the reference ",
      "implementation."
    )
  )
})
