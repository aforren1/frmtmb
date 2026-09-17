# Run every assertion of brms 2.23.0's tests.families.R, and the family
# and response assertions of tests.brm.R and tests.standata.R this lane
# is about, against ONE frmtmb library, and write a ledger.
#
#   Rscript dev/famlink-port-ledger.R base
#   Rscript dev/famlink-port-ledger.R lane
#
# The assertions are PARSED out of brms's own test source, not retyped,
# so the ledger's count is the file's count. Each is evaluated with
# frmtmb attached and brms NOT attached (brms:::internals still resolve
# through the namespace, which is itself a finding: those cannot
# transfer). Output: dev/famlink-port-ledger-<arm>.tsv.
arm <- commandArgs(trailingOnly = TRUE)[1]
lib <- switch(arm,
              base = "C:/Users/adf44/source/r/rellib-r3",
              lane = "C:/Users/adf44/source/r/famlink-lib",
              stop("arm must be base or lane"))
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
suite <- "dev/brms-suite/brms/tests/testthat"
stopifnot(dir.exists(suite))

# (file, test_that label) pairs whose assertions are in scope. The whole
# of tests.families.R; one block each of the other two, and within those
# only the assertions named in `keep`.
blocks <- list(
  list(file = "tests.families.R", label = NULL, keep = NULL),
  list(file = "tests.brm.R", label = NULL,
       keep = "supported link for family|not a supported family"),
  list(file = "tests.standata.R",
       label = "standata rejects incorrect response variables",
       keep = "ordered factors")
)

assertions <- list()
for (b in blocks) {
  exprs <- parse(file.path(suite, b$file), keep.source = FALSE)
  for (e in exprs) {
    if (!is.call(e) || !identical(e[[1]], quote(test_that))) next
    label <- paste(deparse(e[[2]]), collapse = "")
    if (!is.null(b$label) && !grepl(b$label, label, fixed = TRUE)) next
    body <- as.list(e[[3]])[-1]
    setup <- list()
    for (st in body) {
      txt <- paste(deparse(st, width.cutoff = 500L), collapse = " ")
      if (is.call(st) && grepl("^expect_", deparse(st[[1]]))) {
        if (!is.null(b$keep) && !grepl(b$keep, txt)) next
        assertions[[length(assertions) + 1L]] <- list(
          file = b$file, block = label, text = txt, expr = st,
          setup = setup)
      } else {
        setup[[length(setup) + 1L]] <- st
      }
    }
  }
}

# brm() and standata() have no frmtmb spelling. These two shims are
# the TRANSLATION, stated here so the ledger says what was run: brm()
# is frm(), and standata() is frm(dry_run = "frame"), the frame
# assembly where brms's standata() raises the same refusals.
shims <- list(
  brm = function(formula, data, family = NULL, ...) {
    frm(formula, data, family = family)
  },
  standata = function(formula, data, family = NULL, ...) {
    frm(formula, data, family = family, dry_run = "frame")
  }
)

run_one <- function(a) {
  env <- list2env(shims, parent = globalenv())
  # a setup line that fails (mixture(nmix =) is not an frmtmb argument)
  # must not take the block's later, independent assertions with it
  for (s in a$setup) try(eval(s, env), silent = TRUE)
  res <- tryCatch({
    suppressMessages(eval(a$expr, env))
    ""
  }, expectation_failure = function(e) conditionMessage(e),
     error = function(e) paste("ERROR:", conditionMessage(e)))
  # An expect_error() passes on ANY error, including the one R raises
  # for a function or a symbol that does not exist, so its pass says
  # nothing until the error it caught is known. Record that error.
  inner <- ""
  if (is.call(a$expr) && identical(a$expr[[1]], quote(expect_error))) {
    inner <- tryCatch({
      suppressMessages(eval(a$expr[[2]], env))
      ""
    }, error = function(e) conditionMessage(e))
  }
  list(res = res, inner = inner)
}

# An error of this kind means the assertion never reached frmtmb code.
vacuous_pattern <- "could not find function|object '[^']*' not found"

out <- data.frame(idx = seq_along(assertions),
                  file = vapply(assertions, `[[`, "", "file"),
                  assertion = vapply(assertions, `[[`, "", "text"),
                  pass = NA, vacuous = NA, detail = "", caught = "",
                  stringsAsFactors = FALSE)
for (i in seq_along(assertions)) {
  m <- run_one(assertions[[i]])
  out$pass[i] <- !nzchar(m$res)
  out$vacuous[i] <- grepl(vacuous_pattern, m$res) ||
    (out$pass[i] && grepl(vacuous_pattern, m$inner))
  out$detail[i] <- substr(gsub("[\r\n\t]+", " ", m$res), 1, 200)
  out$caught[i] <- substr(gsub("[\r\n\t]+", " ", m$inner), 1, 200)
}
write.table(out, paste0("dev/famlink-port-ledger-", arm, ".tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("arm", arm, ":", sum(out$pass & !out$vacuous), "of", nrow(out),
    "brms assertions pass genuinely;", sum(out$pass & out$vacuous),
    "more pass only on a missing function or object\n")
print(table(file = out$file, genuine = out$pass & !out$vacuous))
