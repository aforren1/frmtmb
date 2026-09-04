# The core/extensions boundary. The protocol is complete.
#
# This test began as a growth-only ratchet: dev/structured-family-
# protocol.md replaced the per-family branches in core with one
# `fam$structure` slot, and the pinned inventory measured the distance
# left to go, shrinking release by release. The extraction round
# emptied it. The ODE seam left for extensions/frmtmb.ode, the sampling
# surface for extensions/frmtmb.sample, and step 10 took the two
# structured families themselves to extensions/frmtmb.latent.
#
# So there is nothing left to ratchet. Every policed token now has an
# empty list of allowed homes: no core file may name hmm, lca, their
# old frame slots, or the ODE symbols, anywhere, at all. This is the
# boundary, not a countdown to it, and a hit is a regression rather
# than debt not yet paid off.
#
# `mixture` is deliberately NOT policed. It stays in core as the
# protocol's reference implementation of the structured-family seam.

# The token vocabulary and the scanner, shared by the boundary
# assertion below and by the positive control that proves the scanner
# can still see a hit.
boundary_tokens <- c("hmm", "lca", "mix_g", "hmm_g", "frm_ode", "RTMBode")

# R's own \b counts `_` as a word character, so `\bhmm\b` would miss
# hmm_check_fit() - exactly the symbols this test exists to see. The
# boundary is therefore spelled against letters and digits only, which
# still keeps `lca` out of `allcaps` and `hmm` out of a longer word.
boundary_scan <- function(rdir) {
  files <- sort(basename(list.files(rdir, pattern = "\\.R$")))
  out <- integer(0)
  for (tok in boundary_tokens) {
    pat <- paste0("(?<![A-Za-z0-9])", tok, "(?![A-Za-z0-9])")
    for (f in files) {
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

# positive identification of the package SOURCE tree, not just an R/
# directory: one CI layout offered an existing-but-empty ../../R, and
# the guard must fail closed (skip) rather than open (assert nothing).
# This is also what makes the test a no-op under `R CMD check` on a
# built tarball, where the tests run beside an INSTALLED package and
# there is no R/ to read.
core_source_dir <- function() {
  rdir <- testthat::test_path("..", "..", "R")
  desc <- testthat::test_path("..", "..", "DESCRIPTION")
  ok <- file.exists(desc) &&
    any(trimws(readLines(desc, n = 5L)) == "Package: frmtmb") &&
    dir.exists(rdir) &&
    file.exists(file.path(rdir, "objective.R"))
  if (ok) rdir else NULL
}

test_that("no core file names a structured family or an ODE symbol", {
  rdir <- core_source_dir()
  skip_if_not(!is.null(rdir),
              "package sources are not available (installed-package run)")

  found <- boundary_scan(rdir)
  hits <- sprintf("  %s: %d", names(found), found)
  expect_identical(
    names(found), NULL,
    info = paste0(
      "A core file names a structured family or an ODE symbol:\n",
      paste(hits, collapse = "\n"),
      "\n\nThe boundary is zero, with no exempt file left in core. Put ",
      "the branch behind the family's own fam$structure slot, register ",
      "it through the exported seams (frmtmb_register_frame_check, ",
      "frmtmb_register_compat), or move the code into the extension ",
      "package that owns it: extensions/frmtmb.latent for hmm and lca, ",
      "extensions/frmtmb.ode for the ODE symbols. mixture() is not ",
      "policed - it stays in core as the reference implementation."
    )
  )
})

test_that("the boundary scanner still sees a token where one exists", {
  # A scanner that matched nothing anywhere would pass the assertion
  # above vacuously. Until step 10 the control read the exempt family
  # homes inside core; they are gone, so it reads them where they now
  # live. This is the only reason this file knows the extension exists.
  #
  # The control is what skips when extensions/ is absent - never the
  # boundary assertion, which must fail closed wherever core sources
  # are readable. Absent extensions/ happens two ways: a tarball check,
  # where .Rbuildignore stripped it (and there the source-tree guard
  # has already skipped everything anyway), and a sparse or partial
  # checkout of the monorepo.
  skip_if_not(!is.null(core_source_dir()),
              "package sources are not available (installed-package run)")
  extdir <- testthat::test_path("..", "..", "extensions", "frmtmb.latent",
                                "R")
  skip_if_not(dir.exists(extdir) && file.exists(file.path(extdir, "hmm.R")),
              "extensions/frmtmb.latent is not in this tree")

  control <- boundary_scan(extdir)
  expect_gt(sum(control), 0L)
  # and specifically for the two tokens that used to be exempt in core,
  # so the control cannot be satisfied by an unrelated match
  expect_true(any(grepl("\\|hmm$", names(control))))
  expect_true(any(grepl("\\|lca$", names(control))))
})
