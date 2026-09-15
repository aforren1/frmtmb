# genrev round 2: is an incomplete `posterior` NO WORSE on this build
# than on base, in a construction other than the lane's?
#
# The lane's construction exports only `ndraws`, which breaks BASE too.
# This one exports every posterior generic the BASE NAMESPACE registers
# on, as real generics, and omits only `as_draws_rvars`: the shape of a
# posterior older than its rvars release.  The lane added ten new
# S3method(posterior::..., frmtmb_fit) directives, and one of them is
# on `as_draws_rvars`.
#
#   Rscript genrev-r2-oldposterior.R <LIB> <order: first|after>
av <- commandArgs(trailingOnly = TRUE)
LIB <- av[1]; order <- av[2]
exps <- c("as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
          "as_draws_matrix", "nchains", "ndraws", "niterations",
          "nvariables", "variables")
lib <- file.path(tempdir(), "oldposterior")
unlink(lib, recursive = TRUE)
d <- file.path(lib, "posterior", "R")
dir.create(d, recursive = TRUE, showWarnings = FALSE)
writeLines(c("Package: posterior", "Version: 0.9.0",
             "Title: Older posterior shape", "Author: t",
             "Maintainer: t <t@t.tt>", "Description: t.",
             "License: GPL-2", "Built: R 4.6.1; ; ; windows"),
           file.path(lib, "posterior", "DESCRIPTION"))
writeLines(sprintf("export(%s)", paste(exps, collapse = ",")),
           file.path(lib, "posterior", "NAMESPACE"))
writeLines(sprintf("%s <- function(x, ...) UseMethod('%s')", exps, exps),
           file.path(d, "posterior"))
.libPaths(c(lib, LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
say <- function(...) cat(sprintf(...))
say("LIB %s  order %s\n", LIB, order)
step <- function(label, expr) {
  r <- tryCatch({ suppressMessages(expr); "ok" },
                error = function(e) paste("FAILED:", substr(conditionMessage(e), 1, 90)))
  say("  %-34s %s\n", label, r)
}
if (order == "first") {
  step("loadNamespace(old posterior)", loadNamespace("posterior"))
  step("library(frmtmb)", library(frmtmb))
} else {
  step("library(frmtmb)", library(frmtmb))
  step("loadNamespace(old posterior)", loadNamespace("posterior"))
}
say("  posterior version loaded: %s\n",
    if (isNamespaceLoaded("posterior")) as.character(packageVersion("posterior")) else "none")
cat("GENREVDONE\n")
