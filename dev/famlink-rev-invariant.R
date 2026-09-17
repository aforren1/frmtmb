## Reviewer check for lane wt-famlink, priority 1: do brms's derived
## link fields on the family the fit carries agree with the links the
## fit actually applies?
##
## fit_assembled() is traced so that every model a test file fits is
## checked at the point the objective is built: for each linear
## predictor in frame$linpreds (the link objective.R and predict.R call),
## the family in spec$responses is read for `link` (dpar mu) or
## `link_<dpar>`, and the two names are compared. A second comparison
## recomputes the derived fields from `links` and asks whether the
## stored ones are what family_link_fields() would write now.
##
## Usage: Rscript dev/famlink-rev-invariant.R <pkgdir> <filter>
## Appends to dev/famlink-rev-invariant-log.tsv. The counts line proves
## the tracer ran (a check that ran zero times is not a pass).
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; filt <- args[2]
ARM <- "lane"
source("dev/famlink-rev-common.R")
suppressMessages(library(testthat))
root <- normalizePath(".")
logf <- file.path(root, "dev", "famlink-rev-invariant-log.tsv")
counter <- new.env()
counter$fits <- 0L; counter$cmp <- 0L; counter$bad <- 0L

check_spec <- function(spec, frame) {
  counter$fits <- counter$fits + 1L
  rec <- function(resp, dpar, field, used, kind) {
    counter$bad <- counter$bad + 1L
    cat(paste(filt, resp, dpar, kind, field, used, sep = "\t"), "\n",
        file = logf, append = TRUE, sep = "")
  }
  for (rn in names(spec$responses)) {
    fam <- spec$responses[[rn]]$family
    if (!inherits(fam, "frmtmb_family")) next
    u <- unclass(fam)
    fresh <- unclass(frmtmb:::family_link_fields(fam))
    for (nm in union(grep("^link", names(u), value = TRUE),
                     grep("^link", names(fresh), value = TRUE))) {
      if (nm %in% c("links", "linkfun", "linkinv")) next
      counter$cmp <- counter$cmp + 1L
      if (!identical(u[[nm]], fresh[[nm]])) {
        rec(rn, nm, paste(format(u[[nm]]), collapse = ","),
            paste(format(fresh[[nm]]), collapse = ","), "stale_vs_links")
      }
    }
    # linkinv on the family against the mu linkinv the objective uses
    for (lp in frame[["linpreds"]]) {
      if (!identical(lp[["resp"]], rn)) next
      dp <- lp[["dpar"]]
      used <- lp[["link"]][["name"]] %||% NA_character_
      field <- if (identical(dp, "mu")) u[["link"]] else
        u[[paste0("link_", dp)]]
      counter$cmp <- counter$cmp + 1L
      if (is.null(field)) {
        if (!identical(fam$type, "categorical") &&
            !grepl("^mu[0-9]*$", dp)) {
          rec(rn, dp, "<absent>", used, "absent")
        }
        next
      }
      if (identical(dp, "mu") && !is.null(u[["ord_link"]])) next
      if (!identical(field, used)) rec(rn, dp, field, used, "fit_uses")
      if (identical(dp, "mu") && is.function(u[["linkinv"]]) &&
          is.function(lp[["link"]][["linkinv"]])) {
        z <- c(-2, -0.3, 0.1, 1.7)
        a <- suppressWarnings(u[["linkinv"]](z))
        b <- suppressWarnings(lp[["link"]][["linkinv"]](z))
        if (!isTRUE(all.equal(a, b))) rec(rn, dp, "linkinv", used,
                                          "linkinv_differs")
      }
    }
  }
}
assign("check_spec", check_spec, envir = globalenv())
suppressMessages(trace("fit_assembled", where = asNamespace("frmtmb"),
      tracer = quote(tryCatch(check_spec(spec, frame), error = function(e)
        cat("TRACER ERROR", conditionMessage(e), "\n"))),
      print = FALSE))

if (!identical(pkgdir, ".")) {
  suppressMessages(library(basename(pkgdir), character.only = TRUE))
}
this_pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
setwd(file.path(pkgdir, "tests", "testthat"))
res <- testthat::test_dir(".", filter = filt, package = this_pkg,
                          reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
line <- sprintf("REVINV %s %s fits %d comparisons %d mismatches %d | pass %d fail %d error %d skip %d",
                pkgdir, filt, counter$fits, counter$cmp, counter$bad,
                sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped))
cat(line, "\n")
cat(line, "\n", file = file.path(root, "dev", "famlink-rev-invariant-counts.txt"),
    append = TRUE, sep = "")
