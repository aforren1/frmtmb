# `[[` on the hazard containers. Keyed on the container's NAME.
#
# THE RULE. In package sources, a `$` is a hit when the token
# immediately to its left is a SYMBOL whose name is in
# `hazard_containers` below. Nothing else is a hit: `$` is not policed
# in general, only on these names. So `fit$obj`, `x$fit` and
# `object$draws` stay exactly as they are, and a chain is judged one
# link at a time by its own left-hand symbol, which is why
# `fit$frame[["par_template"]]` is right and `fit$frame$par_template`
# is not.
#
# WHY. `$` on a list falls back to PARTIAL matching when no name
# matches exactly, so a slot that is absent does not read `NULL`, it
# reads a neighbor whose name it is a prefix of. The containers below
# are keyed by an open or prefix-colliding vocabulary, and the package
# has been bitten three times: `pars$b` read `beta`, `ctx$mix` read
# `mix_g`, and at v0.49 `par_template$theta` read `thetaac` on a model
# with residual autocorrelation and no random effects, so get_prior()
# offered rows for parameters the model did not have. `[[` has no
# fallback and returns NULL, which is what every one of those call
# sites already assumed it was getting.
#
# WHAT IS NOT LISTED. Environments are absent on purpose: `$` on an
# environment is already an exact lookup with no partial fallback, so
# frmtmb_aterm_registry, frmtmb_frame_checks and hyp_shadow_state are
# safe as written. `covstruct_registry` is a plain list and is listed.
# Structs with a closed slot vocabulary (`resp`, `rspec`, `spec`) and
# the user-facing S3 objects (`fit`, `x`, `object`) are not listed
# either: their names are fixed, exact and documented.
#
# CONSEQUENCE. These names are reserved. Binding one of them to
# something that is not the container it names is itself a hit, and
# the fix is to rename the local, not to exempt it.
#
# dev/bracket-sweep.md holds the inventory this list was drawn from,
# with the confirmed-versus-class split at the time of the sweep.

hazard_containers <- c(
  # parameter and estimate lists: theta/thetaac/thetar, b/beta/betad
  "pars", "tpl", "par_template", "est", "estimates",
  # aterm containers: cens/cens_y2, mi/mi_sd, se/se_sigma
  "aterms", "av",
  # random-effects blocks: aux_D/aux_D2, aux_Q/aux_Qk
  "blk", "bk",
  # structure blocks, and the linear-predictor, autocor and eta-design
  # blocks: par/param_colnames, smooth/smooths, p/patterns, n/nonest
  "block", "blocks", "lp", "ac", "ed",
  # the frame itself: y/y_levels, linpred/linpreds
  "frame",
  # family objects: family/family_finalize, mix/mix_groups, sim/sim_ctx
  "fam", "family",
  # simulator contexts, and the distributional-parameter containers,
  # both keyed by names the user supplies
  "ctx", "dpars", "dpv", "dp",
  # registries that are lists rather than environments: us/us_t
  "covstruct_registry"
)

# The scanner. Parse data rather than a regexp, because the left-hand
# token is what decides the hit and only the parser can tell a `$` in
# code from one inside a string or a comment.
bracket_scan <- function(rdir) {
  out <- list()
  for (f in sort(list.files(rdir, pattern = "\\.R$", full.names = TRUE))) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    pd <- pd[pd$terminal & pd$token != "COMMENT", , drop = FALSE]
    pd <- pd[order(pd$line1, pd$col1), , drop = FALSE]
    hit <- which(pd$token == "'$'")
    for (i in hit) {
      if (i <= 1L || pd$token[i - 1L] != "SYMBOL") next
      if (!(pd$text[i - 1L] %in% hazard_containers)) next
      rhs <- if (i < nrow(pd)) pd$text[i + 1L] else "?"
      out[[length(out) + 1L]] <- sprintf(
        "  %s:%d: %s$%s", basename(f), pd$line1[i], pd$text[i - 1L], rhs)
    }
  }
  unlist(out, use.names = FALSE)
}

bracket_advice <- paste0(
  "\n\nUse `[[\"name\"]]`. `$` partial-matches on a list, and every ",
  "container above has names that are prefixes of other names it can ",
  "hold, so the read succeeds with the wrong slot instead of ",
  "returning NULL. If the symbol is not one of these containers, ",
  "rename the local: the names are reserved. See the header of this ",
  "file for the rule and dev/bracket-sweep.md for the inventory."
)

# positive identification of the package SOURCE tree, not just an R/
# directory: one CI layout offered an existing-but-empty ../../R, and
# the guard must fail closed (skip) rather than open (assert nothing).
# This is also what makes the test a no-op under `R CMD check` on a
# built tarball, where the tests run beside an INSTALLED package and
# there is no R/ to read.
bracket_source_dir <- function() {
  rdir <- testthat::test_path("..", "..", "R")
  desc <- testthat::test_path("..", "..", "DESCRIPTION")
  ok <- file.exists(desc) &&
    any(trimws(readLines(desc, n = 5L)) == "Package: frmtmb") &&
    dir.exists(rdir) &&
    file.exists(file.path(rdir, "objective.R"))
  if (ok) rdir else NULL
}

test_that("no core file reads a hazard container with `$`", {
  rdir <- bracket_source_dir()
  skip_if_not(!is.null(rdir),
              "package sources are not available (installed-package run)")

  found <- bracket_scan(rdir)
  expect_identical(
    found, NULL,
    info = paste0("A core file reads a hazard container with `$`:\n",
                  paste(found, collapse = "\n"), bracket_advice)
  )
})

test_that("no extension file reads a hazard container with `$`", {
  # The extensions are policed from here rather than from a copy of
  # this file in each of them, so the container list has one home. The
  # boundary test already reaches across the monorepo the same way.
  # The price is that the guard runs in core's suite and not in each
  # extension's own `R CMD check`; a cloned copy would trade a single
  # container list for five that drift, which is the worse bargain.
  #
  # This is the assertion that skips when extensions/ is absent - never
  # the core one above, which must fail closed wherever core sources
  # are readable. Absent extensions/ happens two ways: a tarball check,
  # where .Rbuildignore stripped it (and there the source-tree guard
  # has already skipped everything anyway), and a sparse or partial
  # checkout of the monorepo.
  skip_if_not(!is.null(bracket_source_dir()),
              "package sources are not available (installed-package run)")
  extroot <- testthat::test_path("..", "..", "extensions")
  dirs <- if (dir.exists(extroot))
    list.dirs(extroot, recursive = FALSE) else character(0)
  dirs <- dirs[dir.exists(file.path(dirs, "R"))]
  skip_if_not(length(dirs) > 0L, "no extensions in this tree")

  for (d in dirs) {
    found <- bracket_scan(file.path(d, "R"))
    expect_identical(
      found, NULL,
      info = paste0("An ", basename(d),
                    " file reads a hazard container with `$`:\n",
                    paste(found, collapse = "\n"), bracket_advice)
    )
  }
})

test_that("the bracket scanner sees a hit where one exists", {
  # The two assertions above are satisfied by an empty tree, so the
  # scanner has to be shown a hit. Unlike the boundary test there is no
  # corner of the monorepo left that still has one - the sweep took
  # them all - so the control writes its own.
  dir <- withr::local_tempdir()
  writeLines(c(
    "f <- function(frame, fit, resp) {",
    "  a <- fit$obj                  # legal: not a hazard container",
    "  b <- resp$resp_name           # legal: closed struct",
    "  d <- \"frame$y\"                # a string, not code",
    "  # frame$y                     # a comment, not code",
    "  e <- frame[[\"y\"]]              # already correct",
    "  c(a, b, d, e, frame$y_levels)",
    "}"), file.path(dir, "control.R"))

  found <- bracket_scan(dir)
  expect_length(found, 1L)
  expect_match(found, "control\\.R:7: frame\\$y_levels", fixed = FALSE)
})
