# Lane optrap: the ported brms suite's six fixtures (helper-brms-suite.R),
# each fitted once: objective, code, message, units.
#   Rscript dev/optrap-fixtures.R <lib or ""> [k ...]
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib", dirname(find.package("frmtmb")), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
src <- readLines("tests/testthat/helper-brms-suite.R")
i <- grep("^brms_fixture_data <- function", src)
j <- grep("^# One fit per fixture per process", src)
eval(parse(text = src[i:(j - 1L)]), envir = env)
ks <- if (length(a) > 1) as.integer(a[-1]) else 1:6
for (k in ks) {
  spec <- env$brms_fixture_spec(k)
  d <- env$brms_fixture_data(k)
  if (!is.null(spec$data)) d <- spec$data(d)
  args <- list(formula = spec$formula, data = d)
  if (!is.null(spec$family)) args$family <- spec$family
  w <- character()
  f <- tryCatch(withCallingHandlers(do.call(frmtmb::frm, args),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("fixture %d: ERROR %s\n", k, conditionMessage(f)))
    next
  }
  cat(sprintf("fixture %d: code %d (%s) objective %.10f | units %s\n", k,
              f$opt$convergence, f$opt$message, f$opt$objective,
              paste(unique(signif(f$par_units %||% 1, 3)), collapse = ",")))
  for (x in w) cat("   W:", substr(x, 1, 140), "\n")
}
