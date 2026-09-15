source("C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-sub.R")

# the detector was too strict: `{ UseMethod("x") }` is a braced body
# with ONE statement and carries no work. Strip one brace layer and
# require exactly one statement, which must be the UseMethod call.
sub1("dev/generics-bodies.R",
paste0("  d <- deparse(body(f))\n",
       "  if (length(d) != 1L || !grepl(\"^UseMethod\\\\(\", trimws(d))) {\n",
       "    bad <- c(bad, sprintf(\"%s: %s\", g,\n",
       "                          paste(trimws(d), collapse = \" \")))\n",
       "  }\n"),
paste0("  b <- body(f)\n",
       "  # a single `{ UseMethod(...) }` carries no work; more than one\n",
       "  # statement inside the braces does\n",
       "  if (is.call(b) && identical(b[[1L]], as.name(\"{\"))) {\n",
       "    b <- if (length(b) == 2L) b[[2L]] else b\n",
       "  }\n",
       "  ok <- is.call(b) && identical(b[[1L]], as.name(\"UseMethod\"))\n",
       "  if (!ok) {\n",
       "    bad <- c(bad, sprintf(\"%s: %s\", g,\n",
       "                          paste(trimws(deparse(body(f))),\n",
       "                                collapse = \" \")))\n",
       "  }\n"))

# and the same guard, in the suite, where it will keep firing
p <- "tests/testthat/test-generic-collision.R"
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
blk <- paste0(
'\n',
'test_that("no shared generic carries work in its own body", {\n',
'  # The structural guard for a defect that cost 8 assertions and was\n',
'  # invisible to a one-file-per-process run.\n',
'  #\n',
'  # frmtmb\'s exported generic for any of these names is REPLACEABLE:\n',
'  # it resolves to the owner\'s generic whenever the owner is loaded,\n',
'  # and the owner\'s generic is a bare UseMethod(). So anything frmtmb\n',
'  # puts in its own generic runs only when no owner is there.\n',
'  # `hypothesis()` used to arm the reserved-name shadowing note in\n',
'  # the generic, and under R CMD check, which runs the whole suite in\n',
'  # ONE process where an earlier file had loaded brms, the note\n',
'  # stopped firing and test-naming-collisions.R lost 8 assertions.\n',
'  # Work of that kind belongs in the METHOD.\n',
'  #\n',
'  # Run in a child with NO owner loaded, so every binding hands back\n',
'  # its own fallback rather than the owner\'s generic.\n',
'  skip_on_cran()\n',
'  out <- run_child(c(\n',
'    "suppressMessages(library(frmtmb))",\n',
'    "ns <- asNamespace(\'frmtmb\')",\n',
'    "tab <- get(\'frm_generic_owners\', envir = ns)",\n',
'    "loaded <- intersect(unique(unlist(tab)), loadedNamespaces())",\n',
'    "cat(\'OWNERSLOADED:\', length(loaded), \'\\\\n\')",\n',
'    "bad <- character()",\n',
'    "for (g in names(tab)) {",\n',
'    "  f <- tryCatch(get(g), error = function(e) NULL)",\n',
'    "  if (!is.function(f)) { bad <- c(bad, g); next }",\n',
'    "  b <- body(f)",\n',
'    "  if (is.call(b) && identical(b[[1L]], as.name(\'{\')))",\n',
'    "    b <- if (length(b) == 2L) b[[2L]] else b",\n',
'    "  if (!(is.call(b) && identical(b[[1L]], as.name(\'UseMethod\'))))",\n',
'    "    bad <- c(bad, g)",\n',
'    "}",\n',
'    "cat(\'NAMES:\', length(tab), \'\\\\n\')",\n',
'    "cat(\'WITHWORK:\', paste(bad, collapse = \',\'), \'\\\\n\')",\n',
'    "cat(\'CHILDOK\\\\n\')"))\n',
'  expect_match(out, "CHILDOK", fixed = TRUE)\n',
'  # the guard is only a guard if the fallbacks are what was inspected\n',
'  expect_equal(parse_field(out, "OWNERSLOADED"), "0")\n',
'  expect_gt(as.integer(parse_field(out, "NAMES")), 20L)\n',
'  expect_equal(parse_field(out, "WITHWORK"), character())\n',
'})\n')
f <- file.path(root, p)
txt <- rawToChar(readBin(f, "raw", file.info(f)$size))
writeBin(charToRaw(paste0(txt, blk)), f)
cat("appended the bare-generic guard\n")
cat("DONE\n")
