# Reviewer, lane wt-conditions: can the core census pass while a bare
# condition call exists? Each case copies R/, inst/, DESCRIPTION and the
# census test into a scratch tree, plants one site, and runs the census.
# A case that PASSES with a plant is a miss.
#   Rscript dev/conditions-rev-census.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
scratch <- file.path(root, "dev/conditions-rev-census")
unlink(scratch, recursive = TRUE)
dir.create(scratch)

make_tree <- function(name) {
  d <- file.path(scratch, name)
  dir.create(file.path(d, "tests", "testthat"), recursive = TRUE)
  file.copy(file.path(root, "R"), d, recursive = TRUE)
  file.copy(file.path(root, "inst"), d, recursive = TRUE)
  file.copy(file.path(root, "DESCRIPTION"), d)
  file.copy(file.path(root, "tests/testthat/test-conditions-census.R"),
            file.path(d, "tests", "testthat"))
  d
}
append_to <- function(path, lines) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  cat(lines, file = path, sep = "\n", append = TRUE)
}

plants <- list(
  control_unplanted = function(d) NULL,
  control_bare_in_R = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() stop("planted")'),
  lowercase_r = function(d)
    append_to(file.path(d, "R/zz-plant.r"),
              'zz_plant <- function() stop("planted")'),
  # R CMD INSTALL collates R/windows/*.R on Windows and R/unix/*.R
  # elsewhere, and accepts .S, .q and .s as R code
  subdir_windows = function(d)
    append_to(file.path(d, "R/windows/zz-plant.R"),
              'zz_plant <- function() stop("planted")'),
  subdir_unix = function(d)
    append_to(file.path(d, "R/unix/zz-plant.R"),
              'zz_plant <- function() stop("planted")'),
  ext_S = function(d)
    append_to(file.path(d, "R/zz-plant.S"),
              'zz_plant <- function() stop("planted")'),
  ext_q = function(d)
    append_to(file.path(d, "R/zz-plant.q"),
              'zz_plant <- function() stop("planted")'),
  sysdata = function(d) {
    zz_plant <- function() stop("planted")
    environment(zz_plant) <- globalenv()
    save(zz_plant, file = file.path(d, "R/sysdata.rda"))
  },
  # a second site under an exemption's file | owner | call key
  exemption_broad = function(d) {
    f <- file.path(d, "R/fit.R")
    x <- readLines(f)
    # optimizer_from_best()'s exempt `    stop(e)` line, duplicated
    i <- grep("^    stop[(]e[)]$", x)[1L]
    x <- append(x, "    stop(e)", after = i - 1L)
    writeLines(x, f)
  },
  exemption_broad_frm_stop = function(d) {
    f <- file.path(d, "R/conditions.R")
    x <- readLines(f)
    # a second stop(cnd) at the top of frm_stop()'s body
    i <- grep("^  # stop[(][)] records the call", x)[1L]
    x <- append(x, "  if (FALSE) stop(cnd)", after = i - 1L)
    writeLines(x, f)
  },
  match_fun = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() match.fun("stop")("planted")'),
  get_string = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() get("warning")("planted")'),
  docall_named_late = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() do.call(args = list("p"), what = "stop")'),
  parse_text = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() eval(parse(text = "stop(\'planted\')"))'),
  str2lang = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() eval(str2lang("warning(\'planted\')"))'),
  deprecated = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() .Deprecated("new_fn")'),
  signal_error = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() base::stop(errorCondition("planted"))'),
  frm_as_value_lapply = function(d)
    append_to(file.path(d, "R/zz-plant.R"),
              'zz_plant <- function() lapply("planted", frm_stop)'),
  inst_subdir_uppercase_S = function(d)
    append_to(file.path(d, "inst/bcm/zz-plant.S"),
              'zz_plant <- function() stop("planted")'),
  vignette_code = function(d)
    append_to(file.path(d, "vignettes/zz-plant.Rmd"),
              c("```{r}", 'stop("planted")', "```"))
)

for (nm in names(plants)) {
  d <- make_tree(nm)
  plants[[nm]](d)
  res <- test_file(file.path(d, "tests/testthat/test-conditions-census.R"),
                   reporter = "silent", package = "frmtmb")
  r <- as.data.frame(res)
  status <- if (sum(r$failed) + sum(r$error) > 0) "FAILS (caught)" else
    if (sum(r$skipped) > 0) "SKIPS" else "PASSES"
  cat(sprintf("%-26s pass %d fail %d error %d skip %d  -> %s\n", nm,
              sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped),
              status))
}
unlink(scratch, recursive = TRUE)
