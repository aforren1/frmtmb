source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")
p <- "tests/testthat/test-generic-collision.R"
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
f <- file.path(root, p)
txt <- rawToChar(readBin(f, "raw", file.info(f)$size))
i <- regexpr('test_that("an owner that exports a NON-generic', txt,
             fixed = TRUE)
if (i < 0) stop("block not found")
txt <- substr(txt, 1, i - 1)

# The first spelling shadowed `posterior`, and that construction could
# not answer the question: frmtmb registers
# `S3method(posterior::as_draws, ...)`, so R's own registerS3methods()
# resolves that name against the shadow at LOAD and frmtmb fails to
# load at all. Measured, and it is PRE-EXISTING: dev/generics-shadowload.R
# gives the identical failure on the base build. `bayesplot` is the
# owner frmtmb registers exactly one method on, so a shadow that
# exports `pp_check` as a plain function satisfies the directive and
# leaves the real question, does frmtmb adopt a NON-generic, standing.
blk <- paste0(
'test_that("an owner that exports a NON-generic does not take the name", {\n',
'  # frm_owner_generic() checks that the owner\'s export is a GENERIC,\n',
'  # not merely a function. Built rather than assumed: a shadow\n',
'  # `bayesplot` whose `pp_check` is a plain function must NOT displace\n',
'  # frmtmb\'s, or frmtmb\'s own method for its own class becomes\n',
'  # unreachable with nothing said.\n',
'  #\n',
'  # bayesplot is the owner frmtmb registers exactly ONE delayed method\n',
'  # on, which is what makes it usable here. Shadowing `posterior` the\n',
'  # same way does not test this at all: R resolves\n',
'  # `S3method(posterior::as_draws, ...)` against the shadow while\n',
'  # LOADING frmtmb, so frmtmb fails to load before any of this runs.\n',
'  # That failure is pre-existing and identical on the base build\n',
'  # (dev/generics-shadowload.R).\n',
'  skip_on_cran()\n',
'  lib <- file.path(tempdir(), "frmtmb-nongeneric-owner")\n',
'  unlink(lib, recursive = TRUE)\n',
'  dir.create(file.path(lib, "bayesplot", "R"), recursive = TRUE,\n',
'             showWarnings = FALSE)\n',
'  on.exit(unlink(lib, recursive = TRUE), add = TRUE)\n',
'  d <- file.path(lib, "bayesplot")\n',
'  writeLines(c("Package: bayesplot", "Version: 99.0",\n',
'               "Title: A deliberately non-generic pp_check",\n',
'               "Author: t", "Maintainer: t <t@t.tt>",\n',
'               "Description: t.", "License: GPL-2",\n',
'               "Built: R 4.6.1; ; ; windows"),\n',
'             file.path(d, "DESCRIPTION"))\n',
'  writeLines("export(pp_check)", file.path(d, "NAMESPACE"))\n',
'  writeLines("pp_check <- function(object, ...) 42",\n',
'             file.path(d, "R", "bayesplot"))\n',
'  out <- run_child(c(\n',
'    sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),\n',
'    "ok <- tryCatch({ loadNamespace(\'bayesplot\'); TRUE },",\n',
'    "               error = function(e) FALSE)",\n',
'    "cat(\'SHADOWLOADED:\', ok, \'\\\\n\')",\n',
'    "cat(\'ISPLAIN:\', !any(grepl(\'UseMethod\',",\n',
'    "      deparse(body(getExportedValue(\'bayesplot\', \'pp_check\'))),",\n',
'    "      fixed = TRUE)), \'\\\\n\')",\n',
'    "suppressMessages(library(frmtmb))",\n',
'    "cat(\'PPCHECK:\', environmentName(topenv(environment(",\n',
'    "      get(\'pp_check\')))), \'\\\\n\')",\n',
'    "cat(\'CHILDOK\\\\n\')"))\n',
'  expect_match(out, "CHILDOK", fixed = TRUE)\n',
'  # the guard is only a guard if the shadow really is there and really\n',
'  # is not a generic\n',
'  expect_equal(parse_field(out, "SHADOWLOADED"), "TRUE")\n',
'  expect_equal(parse_field(out, "ISPLAIN"), "TRUE")\n',
'  expect_equal(parse_field(out, "PPCHECK"), "frmtmb")\n',
'})\n')
writeBin(charToRaw(paste0(txt, blk)), f)
cat("replaced the non-generic block\n")
cat("DONE\n")
