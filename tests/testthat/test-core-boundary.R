# The core/extensions boundary, ratcheted.
#
# dev/structured-family-protocol.md replaces the per-family branches in
# core with one `fam$structure` slot, after which hmm() and lca() can
# live in another package and R/hmm.R, R/lca.R and R/ode.R are the only
# places their names appear. This test measures the distance left to go:
# it counts the family names still reachable from core files and pins
# the count.
#
# The ratchet runs one way. A pinned entry may shrink or disappear -
# that is the refactor working, and the test stays green - but a hit in
# a file/token pair that is not pinned, or a count above what is pinned,
# fails. So a new branch cannot be added to core, and a family name
# cannot spread to a core file it has not already reached.
#
# `mixture` is deliberately NOT policed. It stays in core as the
# protocol's reference implementation of the structured-family seam.
#
# To regenerate the inventory after a refactor step, run the scan below
# over a clean tree and paste the result over `pinned`:
#
#   tokens <- list(hmm = c("hmm.R", "lca.R"), lca = c("hmm.R", "lca.R"),
#                  mix_g = c("hmm.R", "lca.R"), hmm_g = c("hmm.R", "lca.R"),
#                  frm_ode = "ode.R", RTMBode = "ode.R")
#   (the body of boundary_scan() below, with rdir = "R")

test_that("no core file grows a structured-family or ODE reference", {
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

  # each token maps to the files that are allowed to spell it: a
  # family's own home, and R/ode.R for the ODE symbols
  tokens <- list(
    hmm     = c("hmm.R", "lca.R"),
    lca     = c("hmm.R", "lca.R"),
    mix_g   = c("hmm.R", "lca.R"),
    hmm_g   = c("hmm.R", "lca.R"),
    frm_ode = "ode.R",
    RTMBode = "ode.R"
  )

  # R's own \b counts `_` as a word character, so `\bhmm\b` would miss
  # hmm_check_fit() - exactly the branch this test exists to see. The
  # boundary is therefore spelled against letters and digits only, which
  # still keeps `lca` out of `allcaps` and `hmm` out of a longer word.
  boundary_scan <- function(rdir) {
    files <- sort(basename(list.files(rdir, pattern = "\\.R$")))
    out <- integer(0)
    for (tok in names(tokens)) {
      pat <- paste0("(?<![A-Za-z0-9])", tok, "(?![A-Za-z0-9])")
      for (f in setdiff(files, tokens[[tok]])) {
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
    out[order(names(out))]
  }

  # measured after the protocol's steps 1-5 landed (the v0.45.0
  # consolidation re-pinned it from the shrunken v0.44.0 inventory);
  # see the regeneration note above. mix_g is gone entirely, and what
  # remains is the step 6 and 7 residue (hmm_g and the hmm/lca reads
  # in frame, fit, loo, sandwich, predict), compat.R feature-matrix
  # rows, and the two genuine ODE references in interop.R's tmbstan
  # refusal message.
  pinned <- c(
    "compat.R|hmm" = 37L,
    "compat.R|lca" = 39L,
    "fit.R|hmm" = 1L,
    "fit.R|lca" = 2L,
    "frame.R|hmm" = 5L,
    "frame.R|hmm_g" = 4L,
    "frame.R|lca" = 1L,
    "interop.R|RTMBode" = 1L,
    "interop.R|frm_ode" = 1L,
    "loo.R|hmm" = 1L,
    "loo.R|hmm_g" = 1L,
    "methods-draws.R|lca" = 2L,
    "predict.R|lca" = 3L,
    "sandwich.R|hmm" = 2L,
    "sandwich.R|hmm_g" = 1L
  )
  pinned <- pinned[order(names(pinned))]

  found <- boundary_scan(rdir)

  # a scan that matches nothing would pass every assertion below, so
  # check the scanner still sees the seam it is measuring
  expect_gt(sum(found), 0L)

  grew <- character(0)
  for (k in names(found)) {
    was <- if (k %in% names(pinned)) pinned[[k]] else 0L
    if (found[[k]] > was) {
      grew <- c(grew, sprintf("  %s: %d pinned, %d found", k, was,
                              found[[k]]))
    }
  }
  expect_identical(
    grew, character(0),
    info = paste0(
      "A core file gained a structured-family or ODE reference:\n",
      paste(grew, collapse = "\n"),
      "\n\nThis inventory is a one-way ratchet. The structured-family ",
      "protocol (dev/structured-family-protocol.md) shrinks it to ",
      "empty as hmm(), lca() and the ODE seam move behind ",
      "fam$structure; entries may shrink or vanish freely. Nothing may ",
      "grow it: put the new branch behind the family's own slot ",
      "instead, or move the code into R/hmm.R, R/lca.R or R/ode.R. ",
      "mixture() is not policed - it stays in core as the reference ",
      "implementation."
    )
  )
})
