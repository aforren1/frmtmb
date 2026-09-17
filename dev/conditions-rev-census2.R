# Reviewer, lane wt-conditions (round 2): more ways to make the census
# pass with a bare condition call present, on core and on one extension.
# A case that PASSES or SKIPS with a plant is a miss.
#   Rscript dev/conditions-rev-census2.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
scratch <- file.path(root, "dev/conditions-rev-census2")
unlink(scratch, recursive = TRUE)
dir.create(scratch)

make_tree <- function(name, src) {
  d <- file.path(scratch, name)
  dir.create(file.path(d, "tests", "testthat"), recursive = TRUE)
  file.copy(file.path(src, "R"), d, recursive = TRUE)
  if (dir.exists(file.path(src, "inst")))
    file.copy(file.path(src, "inst"), d, recursive = TRUE)
  file.copy(file.path(src, "DESCRIPTION"), d)
  file.copy(file.path(src, "tests/testthat/test-conditions-census.R"),
            file.path(d, "tests", "testthat"))
  d
}
plant <- function(d, code, file = "R/zz-plant.R")
  cat(code, file = file.path(d, file), sep = "\n", append = TRUE)
eam <- file.path(root, "extensions/frmtmb.eam")

cases <- list(
  core_control = list(root, function(d) NULL),
  # the positive identification of a source tree names one file; a
  # rename of that file turns the census into a skip
  core_sentinel_renamed = list(root, function(d) {
    file.rename(file.path(d, "R/objective.R"), file.path(d, "R/objective2.R"))
    plant(d, 'zz <- function() stop("planted")')
  }),
  core_call_string = list(root, function(d)
    plant(d, 'zz <- function() eval(call("stop", "planted"))')),
  core_docall_what_first = list(root, function(d)
    plant(d, 'zz <- function() do.call(what = "stop", list("planted"))')),
  core_value_map = list(root, function(d)
    plant(d, 'zz <- function() Map(warning, "planted")')),
  core_signal_simple_warning = list(root, function(d)
    plant(d, 'zz <- function() .signalSimpleWarning("planted", NULL)')),
  core_exempt_owner_second_site = list(root, function(d)
    plant(d, 'fit_error_context <- function(e) stop(e)')),
  eam_control = list(eam, function(d) NULL),
  eam_base_message = list(eam, function(d)
    plant(d, 'zz <- function() base::message("planted")')),
  eam_sentinel_renamed = list(eam, function(d) {
    file.rename(file.path(d, "R/frmtmb.eam-package.R"),
                file.path(d, "R/package.R"))
    plant(d, 'zz <- function() stop("planted")')
  })
)
for (nm in names(cases)) {
  d <- make_tree(nm, cases[[nm]][[1L]])
  cases[[nm]][[2L]](d)
  res <- test_file(file.path(d, "tests/testthat/test-conditions-census.R"),
                   reporter = "silent", package = "frmtmb")
  r <- as.data.frame(res)
  status <- if (sum(r$failed) + sum(r$error) > 0) "FAILS (caught)" else
    if (sum(r$skipped) > 0) "SKIPS" else "PASSES"
  cat(sprintf("%-32s pass %d fail %d error %d skip %d  -> %s\n", nm,
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
              status))
}
unlink(scratch, recursive = TRUE)
