# Reviewer, lane wt-conditions (recheck 1): legitimate code shapes that
# the round's census rules might flag, run through census_sites() taken
# from the generated core census test, plus the real-tree status.
#   Rscript dev/conditions-rev-census-overreach.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
ex <- parse(file.path(root, "tests/testthat/test-conditions-census.R"),
            keep.source = FALSE)
env <- new.env()
for (e in ex) {
  if (is.call(e) && identical(e[[1L]], as.name("<-"))) eval(e, env)
}
shapes <- c(
  local_var_message = "f <- function(x) { message <- paste('a', x); frm_stop(message) }",
  list_field_message = "f <- function(message) list(message = message)",
  local_get_fn = "f <- function(cache) { get <- function(k) cache[[k]]; get('message') }",
  get_column_env = "f <- function(e) get('warning', envir = e)",
  identical_as_name = "f <- function(h) identical(h, as.name('stop'))",
  call_inspection = "f <- function(e) as.character(e[[1]]) %in% c('stop', 'warning')",
  switch_labels = "f <- function(k) switch(k, message = 1, warning = 2, stop = 3)",
  local_fn_named_call = "f <- function() { call <- function(x) x; call('stop') }",
  parse_formula_text = "f <- function() str2lang('y ~ message(x)')",
  mget_names = "f <- function(e) mget(c('a', 'message'), envir = e)",
  getOption_warn = "f <- function() getOption('warning.length')",
  formals_default_msg = "f <- function(on_fail = c('stop', 'warning')) match.arg(on_fail)",
  R6_like = "f <- function(obj) obj$get('message')",
  dollar_call = "f <- function(x) x$call('stop')")
for (nm in names(shapes)) {
  got <- env$census_sites(parse(text = shapes[[nm]], keep.source = FALSE),
                          "R/x.R")
  cat(sprintf("%-22s %s\n", nm, if (length(got)) paste("FLAGGED:",
      paste(sub("^R/x[.]R [|] f [|] ", "", got), collapse = "; ")) else "clean"))
}
