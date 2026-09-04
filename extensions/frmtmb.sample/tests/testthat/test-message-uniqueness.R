# The property frmtmb asserts on its own R/, asserted here on this
# package's: every condition message TEMPLATE raised in R/ is unique,
# for stop(), warning() and message() alike. A template is the
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
# The property is per PACKAGE, not across the pair. A template here may
# repeat one of frmtmb's without ambiguity, because the two resolve to
# different source trees; what must not happen is two lines of THIS
# package rendering the same message.

test_that("every condition message template in this R/ is unique", {
  # positive identification of the package SOURCE tree, not just an R/
  # directory: one CI layout offered an existing-but-empty ../../R, and
  # the guard must fail closed (skip) rather than open (assert nothing)
  rdir <- testthat::test_path("..", "..", "R")
  desc <- testthat::test_path("..", "..", "DESCRIPTION")
  is_src <- file.exists(desc) &&
    any(trimws(readLines(desc, n = 5L)) == "Package: frmtmb.sample") &&
    dir.exists(rdir) &&
    file.exists(file.path(rdir, "sample.R"))
  skip_if_not(is_src,
              "package sources are not available (installed-package run)")

  collect <- function(kind) {
    out <- character(0)
    walk <- function(e) {
      if (!is.call(e)) return(invisible(NULL))
      if (is.name(e[[1L]]) && identical(as.character(e[[1L]]), kind)) {
        lits <- character(0)
        for (i in seq_along(e)[-1L]) {
          # an empty argument (`x[, 1]`) errors when TOUCHED, not when
          # extracted, so the inspection itself is guarded
          v <- tryCatch(if (is.character(e[[i]])) e[[i]],
                        error = function(err) NULL)
          if (is.character(v)) lits <- c(lits, paste(v, collapse = ""))
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

  for (kind in c("stop", "warning", "message")) {
    msgs <- collect(kind)
    expect_gt(length(msgs), 0, label = paste0("templates found for ", kind))
    dup <- unique(msgs[duplicated(msgs)])
    expect_length(dup, 0)
    if (length(dup)) {
      fail(paste0("duplicated ", kind, "() templates: ",
                  paste(substr(dup, 1, 60), collapse = " | ")))
    }
  }
})
