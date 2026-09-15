# Does an INCOMPLETE owner on the library path stop frmtmb loading,
# and is that this lane's doing?
#
# Found by a test of mine that shadowed `posterior` with a package
# exporting only `ndraws`.  frmtmb then failed to load at all:
#
#   Error: package or namespace load failed for 'frmtmb' in
#   get(genname, envir = envir): object 'as_draws' not found
#
# That is R's own `registerS3methods()` resolving
# `S3method(posterior::as_draws, ...)` against a namespace that does
# not export the name.  The question this script answers is whether
# the BASE build does the same, because `S3method(posterior::as_draws,
# frmtmb_multiple)` is in the base commit's NAMESPACE too.
#
#   Rscript dev/generics-shadowload.R <LIB> <owner> <exports...>
av <- commandArgs(trailingOnly = TRUE)
LIB <- av[1]
own <- av[2]
exps <- av[-(1:2)]
lib <- file.path(tempdir(), paste0("shadow-", own))
unlink(lib, recursive = TRUE)
d <- file.path(lib, own, "R")
dir.create(d, recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste0("Package: ", own), "Version: 99.0",
             "Title: Deliberately incomplete", "Author: t",
             "Maintainer: t <t@t.tt>", "Description: t.",
             "License: GPL-2", "Built: R 4.6.1; ; ; windows"),
           file.path(lib, own, "DESCRIPTION"))
writeLines(sprintf("export(%s)", paste(exps, collapse = ",")),
           file.path(lib, own, "NAMESPACE"))
writeLines(vapply(exps, function(e)
  sprintf("%s <- function(x, ...) 42", e), ""),
  file.path(d, own))

f <- tempfile(fileext = ".R")
writeLines(c(
  sprintf(".libPaths(c(%s, %s, %s, %s, %s))", deparse(lib), deparse(LIB),
          deparse("C:/Users/adf44/source/r/rellib-r3"),
          deparse("C:/Users/adf44/source/r/pinlib"),
          deparse("C:/Users/adf44/AppData/Local/R/win-library/4.6")),
  sprintf("cat('SHADOW %s exports %s\\n')", own,
          paste(exps, collapse = ",")),
  sprintf("ok <- tryCatch({ loadNamespace('%s'); TRUE },", own),
  "               error = function(e) FALSE)",
  "cat('SHADOW_LOADS:', ok, '\\n')",
  "r <- tryCatch({ suppressMessages(library(frmtmb)); 'loaded' },",
  "              error = function(e) paste('FAILED:',",
  "                substr(conditionMessage(e), 1, 70)))",
  "cat('FRMTMB:', r, '\\n')",
  "if (identical(r, 'loaded')) {",
  sprintf("  cat('%s RESOLVES TO:',", exps[1]),
  sprintf("      environmentName(topenv(environment(get('%s')))), '\\n')",
          exps[1]),
  "  set.seed(1); dd <- data.frame(x = rnorm(40))",
  "  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)",
  "  fit <- frm(bf(y ~ x) + gaussian(), data = dd)",
  sprintf("  v <- tryCatch(%s(fit), error = function(e)", exps[1]),
  "                paste('ERR', substr(conditionMessage(e), 1, 40)))",
  sprintf("  cat('%s(fit):', class(v)[1], '\\n')", exps[1]),
  "}",
  "cat('PROBEDONE\\n')"), f)
o <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                              c("--vanilla", shQuote(f)),
                              stdout = TRUE, stderr = TRUE))
cat(paste(grep("^SHADOW|^FRMTMB|RESOLVES TO|\\(fit\\):|^PROBEDONE", o,
               value = TRUE), collapse = "\n"), "\n")
