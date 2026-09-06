# The property frmtmb asserts on its own R/, asserted here on this
# package's: every condition message TEMPLATE raised in R/ is unique,
# for stop(), warning() and message() alike, and for the
# errorCondition() and warningCondition() constructors, which raise
# through stop() and warning() with the message one call further in. A template is the
# concatenation of a call's literal string fragments; runtime
# interpolation (argument names, family names) is not part of it, so
# what this guarantees is that a reported message resolves to one line
# of source, not that two calls can never render the same final text.
# The shared validation helpers are the known case: one template each,
# the argument name filled in at run time.
#
# Counting was retired deliberately: a number in prose is stale the
# commit after it is written, and this test fails at the moment of
# drift instead.
#
# The property is per PACKAGE, not across the repository. A template
# here may repeat one of frmtmb's without ambiguity, because the two
# resolve to different source trees; what must not happen is two lines
# of THIS package rendering the same message.

test_that("every condition message template in this R/ is unique", {
  # positive identification of the package SOURCE tree, not just an R/
  # directory: one CI layout offered an existing-but-empty ../../R, and
  # the guard must fail closed (skip) rather than open (assert nothing)
  rdir <- testthat::test_path("..", "..", "R")
  desc <- testthat::test_path("..", "..", "DESCRIPTION")
  is_src <- file.exists(desc) &&
    any(trimws(readLines(desc, n = 5L)) == "Package: frmtmb.spline") &&
    dir.exists(rdir) &&
    file.exists(file.path(rdir, "curve-cov.R"))
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

  found <- 0L
  pools <- list(error = c("stop", "errorCondition"),
                warning = c("warning", "warningCondition"),
                message = "message")
  for (kind in names(pools)) {
    msgs <- collect(pools[[kind]])
    found <- found + length(msgs)
    dup <- unique(msgs[duplicated(msgs)])
    expect_length(dup, 0)
    if (length(dup)) {
      fail(paste0("duplicated ", kind, " templates: ",
                  paste(substr(dup, 1, 60), collapse = " | ")))
    }
  }
  # the parser's own sanity check, summed over the three pools rather
  # than asserted per pool: this package raises no message() at all, so
  # a per-pool assertion would fail on a true statement about it. Most
  # of what it raises is a refusal, because it names a thing the
  # construction cannot do rather than a thing that went wrong; the
  # warnings are the cases where the answer is still an answer and the
  # reader needs to know how it was reached. The sum guards against a
  # collect() that silently matches nothing, which is what it is for.
  expect_gt(found, 0)
})
