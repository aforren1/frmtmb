# The G5.2a property, asserted rather than counted: every condition
# message TEMPLATE raised in R/ is unique, for stop(), warning() and
# message() alike, and for the errorCondition() and warningCondition()
# constructors, which raise through stop() and warning() with the
# message one call further in. A constructor is pooled with the plain
# form it is raised through, so a warningCondition() template that
# repeated a warning() template would fail here. A template is the concatenation of a call's literal
# string fragments; runtime interpolation (argument names, family
# names) is not part of it, so what this guarantees is that a reported
# message resolves to one line of source, not that two calls can never
# render the same final text. The shared validation helpers are the
# known case: one template each, the argument name filled in at run
# time.
#
# Counting was retired deliberately: a number in prose is stale the
# commit after it is written, and this test fails at the moment of
# drift instead.

test_that("every condition message template in R/ is unique", {
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

  collect <- function(kinds) {
    out <- character(0)
    walk <- function(e) {
      if (!is.call(e)) return(invisible(NULL))
      if (is.name(e[[1L]]) && as.character(e[[1L]]) %in% kinds) {
        lits <- character(0)
        nms <- names(e)
        if (is.null(nms)) nms <- rep("", length(e))
        for (i in seq_along(e)[-1L]) {
          # `class =` on a condition constructor is a character literal
          # and is not message text; neither is `call =`. Skipping them
          # by NAME is what lets warningCondition() into the pool.
          if (nzchar(nms[i]) &&
              nms[i] %in% c("class", "call", "call.", "domain",
                            "appendLF")) {
            next
          }
          # an empty argument (`x[, 1]`) errors when TOUCHED, not when
          # extracted, so the inspection itself is guarded
          v <- tryCatch(if (is.character(e[[i]])) e[[i]],
                        error = function(err) NULL)
          if (is.character(v)) {
            lits <- c(lits, paste(v, collapse = ""))
            next
          }
          # a condition constructor takes ONE message argument, so its
          # fragments are one paste0() deeper than stop()'s are; without
          # this the templates raised through errorCondition() and
          # warningCondition() would be collected as the empty string
          # and drop out of the pool unasserted
          a <- tryCatch(if (is.call(e[[i]])) e[[i]], error = function(err) NULL)
          if (!is.null(a) && is.name(a[[1L]]) &&
              as.character(a[[1L]]) %in% c("paste0", "paste")) {
            for (j in seq_along(a)[-1L]) {
              u <- tryCatch(if (is.character(a[[j]])) a[[j]],
                            error = function(err) NULL)
              if (is.character(u)) lits <- c(lits, paste(u, collapse = ""))
            }
          }
        }
        txt <- paste(lits, collapse = "")
        if (nzchar(txt)) out <<- c(out, txt)
      }
      for (i in seq_along(e)) {
        a <- tryCatch(if (is.call(e[[i]])) e[[i]],
                      error = function(err) NULL)
        if (!is.null(a)) walk(a)
      }
      invisible(NULL)
    }
    for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE)) {
      for (e in parse(f, keep.source = FALSE)) walk(e)
    }
    out
  }

  pools <- list(error = c("stop", "errorCondition"),
                warning = c("warning", "warningCondition"),
                message = "message")
  for (kind in names(pools)) {
    msgs <- collect(pools[[kind]])
    expect_gt(length(msgs), 0, label = paste0("templates found for ", kind))
    dup <- unique(msgs[duplicated(msgs)])
    expect_length(dup, 0)
    if (length(dup)) {
      fail(paste0("duplicated ", kind, " templates: ",
                  paste(substr(dup, 1, 60), collapse = " | ")))
    }
  }
})
