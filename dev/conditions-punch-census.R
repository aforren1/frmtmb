# Lane wt-conditions, punch round 1: the census holes of review MINOR 3,
# planted one at a time into a scratch copy of a tree, with the tree's
# own census test run on each. Every plant must FAIL; the two unplanted
# controls must PASS with nothing skipped. The plants are those of
# dev/conditions-rev-census.R and dev/conditions-rev-census2.R, whose
# logs (census-rerun.txt, census2.txt) record the census before this
# round passing or skipping on each of them.
#   Rscript dev/conditions-punch-census.R > dev/conditions-punch-log/census.txt
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
scratch <- file.path(root, "dev/conditions-punch-census-scratch")
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
plant <- function(d, code, file = "R/zz-plant.R") {
  path <- file.path(d, file)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  cat(code, file = path, sep = "\n", append = TRUE)
}
dup_line <- function(d, file, pattern, line) {
  f <- file.path(d, file)
  x <- readLines(f)
  i <- grep(pattern, x)[1L]
  stopifnot(!is.na(i))
  writeLines(append(x, line, after = i - 1L), f)
}
eam <- file.path(root, "extensions/frmtmb.eam")
bare <- 'zz_plant <- function() stop("planted")'

cases <- list(
  core_control = list(root, function(d) NULL),
  eam_control = list(eam, function(d) NULL),
  # counted exemptions
  core_exempt_second_site_fit = list(root, function(d)
    dup_line(d, "R/fit.R", "^    stop[(]e[)]$", "    stop(e)")),
  core_exempt_second_site_helper = list(root, function(d)
    dup_line(d, "R/conditions.R", "^  stop[(]cnd[)]$",
             "  if (FALSE) stop(cnd)")),
  core_exempt_owner_redefined = list(root, function(d)
    plant(d, 'fit_error_context <- function(e) stop(e)')),
  # files R CMD INSTALL reads
  core_lowercase_r = list(root, function(d)
    plant(d, bare, "R/zz-plant.r")),
  core_subdir_windows = list(root, function(d)
    plant(d, bare, "R/windows/zz-plant.R")),
  core_subdir_unix = list(root, function(d)
    plant(d, bare, "R/unix/zz-plant.R")),
  core_ext_S = list(root, function(d) plant(d, bare, "R/zz-plant.S")),
  core_ext_q = list(root, function(d) plant(d, bare, "R/zz-plant.q")),
  core_ext_lower_s = list(root, function(d)
    plant(d, bare, "R/zz-plant.s")),
  core_inst_subdir_S = list(root, function(d)
    plant(d, bare, "inst/bcm/zz-plant.S")),
  # do.call with `what` named
  core_docall_what_first = list(root, function(d)
    plant(d, 'zz <- function() do.call(what = "stop", list("planted"))')),
  core_docall_what_late = list(root, function(d)
    plant(d, 'zz <- function() do.call(args = list("p"), what = "stop")')),
  # lookups, raisers and built calls
  core_get = list(root, function(d)
    plant(d, 'zz <- function() get("warning")("planted")')),
  core_match_fun = list(root, function(d)
    plant(d, 'zz <- function() match.fun("stop")("planted")')),
  core_call_string = list(root, function(d)
    plant(d, 'zz <- function() eval(call("stop", "planted"))')),
  core_str2lang = list(root, function(d)
    plant(d, 'zz <- function() eval(str2lang("warning(\'planted\')"))')),
  core_parse_text = list(root, function(d)
    plant(d, 'zz <- function() eval(parse(text = "stop(\'planted\')"))')),
  core_deprecated = list(root, function(d)
    plant(d, 'zz <- function() .Deprecated("new_fn")')),
  core_signal_simple_warning = list(root, function(d)
    plant(d, 'zz <- function() .signalSimpleWarning("planted", NULL)')),
  # a renamed sentinel
  core_sentinel_renamed = list(root, function(d) {
    file.rename(file.path(d, "R/objective.R"),
                file.path(d, "R/objective2.R"))
  }),
  eam_sentinel_renamed = list(eam, function(d) {
    file.rename(file.path(d, "R/frmtmb.eam-package.R"),
                file.path(d, "R/package.R"))
  }),
  # the shapes caught before this round, kept as controls on the rule
  core_bare = list(root, function(d) plant(d, bare)),
  eam_base_message = list(eam, function(d)
    plant(d, 'zz <- function() base::message("planted")')),
  core_value_map = list(root, function(d)
    plant(d, 'zz <- function() Map(warning, "planted")'))
)
for (nm in names(cases)) {
  d <- make_tree(nm, cases[[nm]][[1L]])
  cases[[nm]][[2L]](d)
  res <- test_file(file.path(d, "tests/testthat/test-conditions-census.R"),
                   reporter = "silent", package = "frmtmb")
  r <- as.data.frame(res)
  status <- if (sum(r$failed) + sum(r$error) > 0) "FAILS" else
    if (sum(r$skipped) > 0) "SKIPS" else "PASSES"
  want <- if (grepl("_control$", nm)) "PASSES" else "FAILS"
  cat(sprintf("%-34s pass %2d fail %d error %d skip %d  -> %-6s %s\n", nm,
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
              status, if (status == want) "as required" else "WRONG"))
}
unlink(scratch, recursive = TRUE)
