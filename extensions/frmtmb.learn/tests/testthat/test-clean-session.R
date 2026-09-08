## Does `rlddm()` work for a user who loads the package and fits data
## they already have?
##
## IT DID NOT, and nothing in this suite could tell. `rlddm()` requires
## the `dec()` addition term, which reaches frmtmb's registry from
## frmtmb.eam's own `.onLoad()`. A DESCRIPTION `Imports:` entry only
## requires a package to be INSTALLED; it does not load it, and
## `frmtmb.eam::wiener_lpdf()` loads the namespace when the function is
## CALLED, which is after `frm()` has parsed the formula and refused the
## term. Every other rlddm block in this package calls
## `frm_task_simulate()` first, which reaches `ddm_simulate()` and loads
## the namespace as a side effect, so every one of them masked it.
##
## THE TEST THEREFORE RUNS A FRESH R PROCESS. It cannot be written
## in-process: by the time testthat runs, this package is loaded and so
## are its imports, which is precisely the state the bug hides in. The
## file is its own so that no earlier block in the same file can load
## the namespace for it either.

test_that("rlddm() fits from a clean session on a user's own data", {
  skip_on_cran()
  rscript <- file.path(R.home("bin"),
                       if (.Platform$OS.type == "windows") {
                         "Rscript.exe"
                       } else {
                         "Rscript"
                       })
  skip_if_not(file.exists(rscript), "no Rscript to start a clean session")

  # The library paths are written INTO the child script rather than
  # passed through system2(env =), which is documented as unsupported on
  # Windows and there produces no output at all. Without them the child
  # would find whatever is installed system-wide instead of the packages
  # this suite is running against.
  body <- '
suppressPackageStartupMessages(library(frmtmb.learn))
cat("EAM:", "frmtmb.eam" %in% loadedNamespaces(), "\n")
f <- frmtmb::frm_compat_features()
cat("DEC:", "dec" %in% f$key[f$kind == "aterm"], "\n")
set.seed(1)
n <- 120L
d <- data.frame(id = rep(1:6, each = 20), trial = rep(1:20, times = 6),
                rt = runif(n, 0.4, 1.4), choice = rbinom(n, 1, 0.5),
                pay1 = rbinom(n, 1, 0.7), pay2 = rbinom(n, 1, 0.3))
fit <- frmtmb::frm(
  frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
             bs ~ 1, ndt ~ 1, bias ~ 1),
  family = rlddm(subject = id, trial = trial), data = d)
cat("FIT:", is.finite(as.numeric(stats::logLik(fit))), "\n")
'
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(paste0(".libPaths(", paste(deparse(.libPaths()),
                                          collapse = ""), ")"),
               body), f)
  out <- suppressWarnings(system2(rscript, shQuote(f), stdout = TRUE,
                                  stderr = TRUE))
  info <- paste(out, collapse = "\n")

  # the namespace is loaded by the import, before anything is called
  expect_true(any(grepl("^EAM: TRUE", out)), info = info)
  # so the addition term the family requires is in the registry
  expect_true(any(grepl("^DEC: TRUE", out)), info = info)
  # and the fit the user actually wanted goes through. No
  # frm_task_simulate() anywhere above: the data frame is built by hand,
  # which is the one thing every other rlddm test does not do.
  expect_true(any(grepl("^FIT: TRUE", out)), info = info)
  expect_false(any(grepl("dec... is not supported", out)), info = info)
})
