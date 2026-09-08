# `[[` on the hazard containers. Keyed on the container's NAME.
#
# THE RULE, THE LIST AND THE WHY now live in `R/hazard-containers.R`,
# beside the exported scanner `frm_hazard_reads()`, because the same
# rule has to reach seven extensions and a list with two homes has
# none. Read that file first; this one is the SOURCE-TREE half of the
# guard.
#
# In short: a `$` is a hit when the token immediately to its left is a
# SYMBOL whose name is in `frmtmb:::hazard_containers`. Nothing else is
# a hit, so `fit$obj` and `object$draws` stay exactly as they are, and a
# chain is judged one link at a time by its own left-hand symbol. `$` on
# a list partial-matches, so a `$` read of an absent slot on one of
# those containers returns a neighbor instead of NULL.
#
# TWO SCANNERS, ON PURPOSE, and this file pins them against each other.
#
# * `bracket_scan()` below reads the SOURCES, so it reports file and
#   line, it sees top-level code that is not inside a function, and it
#   can be pointed at any directory. It needs the monorepo checkout,
#   and it is a no-op under `R CMD check` on a built tarball, where the
#   tests run beside an installed package and there is no R/ to read.
# * `frm_hazard_reads()` reads the installed NAMESPACE, so it runs
#   wherever the package does, including inside an extension's own
#   `R CMD check`. That is what each extension's own
#   `test-bracket-access.R` calls, in eight lines of assertion,
#   with no copy of the list.
#
# dev/bracket-sweep.md holds the inventory the list was drawn from,
# with the confirmed-versus-class split at the time of the sweep.

hazard_containers <- frmtmb:::hazard_containers

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
  # THE BACKSTOP, not the guard. Each extension now fails its OWN check
  # on a hazard read, through frm_hazard_reads() in its own
  # test-bracket-access.R; frmtmb.coupling passed its check and its
  # suite with 34 of them and was caught here, at the release tally,
  # which is too late. This assertion stays because it reads the
  # SOURCES: it also covers an extension that deleted its own guard
  # file, and top-level code that is not inside a function, which the
  # namespace scan cannot see.
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

# ---------------------------------------------------------------------
# THE NAMESPACE SCANNER, which is the half that runs inside an
# extension's own `R CMD check`. Everything below is about
# frm_hazard_reads(): that it sees a hit, that it refuses the same four
# things bracket_scan() refuses, that the two agree, and that every
# extension in the tree actually calls it.
# ---------------------------------------------------------------------

test_that("frm_hazard_reads() sees a hit, and only a real one", {
  # The same control as the source scanner's, as an environment: a
  # legal container, a closed struct, a `$` inside a string, an
  # already-correct `[[`, and one hit.
  e <- new.env()
  e$f <- function(frame, fit, resp) {
    a <- fit$obj
    b <- resp$resp_name
    d <- "frame$y"
    g <- frame[["y"]]
    c(a, b, d, g, frame$y_levels)
  }
  # a default argument is code that runs in the callee's frame, so it
  # is scanned too
  e$h <- function(pars, k = pars$beta) k
  # a chain is judged one link at a time by its own left-hand name, so
  # the SECOND `$` here is the hit and `fit$frame` is not. This is the
  # idiom the pre-sweep frmtmb.coupling had at coupling.R:185, and a
  # scan that only looked at the outer expression's left side missed it
  # while the source scan caught it.
  e$k <- function(fit) c(fit$frame$data, fit[["frame"]]$par_template)
  # an INNER function's default argument. Formals are a pairlist, not a
  # call, so a walk that descends only into calls steps over it while
  # the source scan reads it; the two scanners have to agree here.
  e$m <- function(x) {
    inner <- function(ctx, j = ctx$mix) j
    inner(x)
  }
  hits <- frmtmb::frm_hazard_reads(e)
  expect_identical(hits, c("f: frame$y_levels", "h: pars$beta",
                           "k: frame$data", "m: ctx$mix"))

  # the reserved-name rule: `family` is a hazard container, and a local
  # bound to something else does not earn an exemption
  e2 <- new.env()
  e2$g <- function(family) family$name
  expect_identical(frmtmb::frm_hazard_reads(e2), "g: family$name")

  expect_error(frmtmb::frm_hazard_reads(c("a", "b")),
               "one package name", fixed = TRUE)
})

test_that("the two scanners agree, and core is clean by both", {
  expect_identical(frmtmb::frm_hazard_reads("frmtmb"), character(0))

  rdir <- bracket_source_dir()
  skip_if_not(!is.null(rdir),
              "package sources are not available (installed-package run)")
  # Both scanners on the same control: same one hit, reported in each
  # one's own vocabulary (file and line for the sources, function name
  # for the namespace). This is what keeps the namespace scan from
  # going blind while the source scan carries the suite.
  dir <- withr::local_tempdir()
  writeLines(c(
    "f <- function(frame, fit) {",
    "  a <- fit$obj",
    "  b <- fit$frame$data",
    "  g <- function(ctx, j = ctx$mix) j",
    "  c(a, b, g(fit), frame[[\"y\"]], frame$y_levels)",
    "}"), file.path(dir, "control.R"))
  e <- new.env()
  source(file.path(dir, "control.R"), local = e)

  src <- bracket_scan(dir)
  ns <- frmtmb::frm_hazard_reads(e)
  expect_length(src, 3L)
  expect_length(ns, 3L)
  # the same three reads, each in its own vocabulary. Measured on the
  # pre-sweep frmtmb.coupling sources, where the source scan reported
  # 34 and this one reported the same 34.
  expect_identical(sub("^.*: ", "", src), sub("^.*: ", "", ns))
  expect_match(src, "frame\\$data", all = FALSE)
  expect_match(ns, "frame\\$data", all = FALSE)
})

test_that("every extension in the tree runs the guard on itself", {
  # The recorded failure this answers: frmtmb.coupling passed its own
  # `R CMD check` and its own suite with 34 hazard reads, because the
  # guard lived only here. A new extension that does not carry the file
  # now fails core's suite, which is the earliest place that can notice
  # it at all.
  skip_if_not(!is.null(bracket_source_dir()),
              "package sources are not available (installed-package run)")
  extroot <- testthat::test_path("..", "..", "extensions")
  dirs <- if (dir.exists(extroot))
    list.dirs(extroot, recursive = FALSE) else character(0)
  dirs <- dirs[dir.exists(file.path(dirs, "R"))]
  skip_if_not(length(dirs) > 0L, "no extensions in this tree")

  for (d in dirs) {
    pkg <- basename(d)
    f <- file.path(d, "tests", "testthat", "test-bracket-access.R")
    expect_true(
      file.exists(f),
      info = paste0(pkg, " has no tests/testthat/test-bracket-access.R, ",
                    "so a hazard read in it would pass its own check. ",
                    "Copy any other extension's: the assertion is ",
                    "eight lines and carries no container list."))
    if (!file.exists(f)) next
    txt <- paste(readLines(f, warn = FALSE), collapse = " ")
    # its own name, not a neighbor's: a copied file that still scans
    # the package it was copied from asserts nothing about this one
    expect_match(txt, paste0("frm_hazard_reads(\"", pkg, "\")"),
                 fixed = TRUE)
  }
})
