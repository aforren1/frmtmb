# One package's whole test suite in ONE R session, the way covr and
# R CMD check run it, in alphabetical or reverse file order, or with
# every extension loaded before the first file.
#
#   Rscript run-session.R <package> [forward|reverse|preload] [test dir]
#
# WHY A SECOND RUNNER BESIDE run-tests.R. That one runs one file per
# process, which isolates each file from what the files before it
# loaded and registered. covr and R CMD check do not: they run every
# file in one session, alphabetically. 0.71.0's test-coverage job went
# red there and nowhere else, because test-families-extreme.R loaded
# frmtmb.eam, whose .onLoad() registers dec(), and test-parse.R assumed
# dec() was unregistered. Reverse order exposes the other direction: a
# file that passes only because no earlier file loaded an extension.
# Preload (alphabetical order, all seven extensions loaded first) is the
# session in which no file loads an frmtmb namespace, so every registry
# difference it reports is the test's own; it is also the session in
# which a test that assumes an extension is absent shows.
#
# Every extension is required to be installed, because a suite that
# cannot load one never meets the registrations it makes.
#
# After each file the runner compares frmtmb's registries, the options,
# the environment variables, the search path and the global environment
# with what they were before the file, and prints a LEAK line for each
# difference. A file must leave every registry as it found it
# (dev/lane-rules.md, "Test isolation").
#
# The caller sets the tier's environment (NOT_CRAN, gates) as for
# run-par.sh. Output, on stdout:
#   lib:     which build of each package the session loaded
#   FILE     a file starting
#   DETAIL   every failure, error, skip and escaped warning, with its test
#   LEAK     a registry or global the file left changed
#   RESULT   per file, in run-tests.R's format
#   SESSION  the totals, last
local({
# The same libraries as run-tests.R: FRMTMB_RELLIB ahead of the release
# library, for a lane's build. Change LIB in both files at a release.
LIB <- "C:/Users/adf44/source/r/rellib-r15"
.libPaths(c(if (nzchar(Sys.getenv("FRMTMB_RELLIB"))) {
  Sys.getenv("FRMTMB_RELLIB")
}, LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))

a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
ord <- if (length(a) >= 2L) a[2] else "forward"
stopifnot(ord %in% c("forward", "reverse", "preload"))
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                 value = TRUE)[1])
root <- normalizePath(file.path(dirname(self), "..", ".."), "/")
dir <- if (length(a) >= 3L) a[3] else if (identical(p, "frmtmb")) {
  file.path(root, "tests", "testthat")
} else {
  file.path(root, "extensions", p, "tests", "testthat")
}
stopifnot(dir.exists(dir))

pkgs <- c("frmtmb", "frmtmb.coupling", "frmtmb.eam", "frmtmb.latent",
          "frmtmb.learn", "frmtmb.ode", "frmtmb.sample", "frmtmb.survival",
          "frmtmb.spatial")
miss <- pkgs[!nzchar(vapply(pkgs, function(q) {
  tryCatch(find.package(q), error = function(e) "")
}, ""))]
if (length(miss)) {
  stop("not installed, and a one-session run needs every extension: ",
       paste(miss, collapse = ", "))
}
for (q in pkgs) {
  cat("lib: ", q, " ", as.character(utils::packageVersion(q)), " ",
      dirname(find.package(q)), "\n", sep = "")
}

suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
if (identical(ord, "preload")) {
  for (q in setdiff(pkgs, "frmtmb")) loadNamespace(q)
}
# the order covr and R CMD check use, from testthat's own file finder
paths <- testthat:::find_test_scripts(dir, full.names = FALSE)
if (identical(ord, "reverse")) paths <- rev(paths)
cat("order: ", ord, ", ", length(paths), " files from ", dir, "\n",
    sep = "")

# What a file can leave behind for the files after it. The frmtmb
# registries are read through the namespace, never copied: a copy would
# be the thing compared against itself. A name the namespace does not
# hold stops the run: `ns$<renamed>` is NULL, and NULL compared with
# NULL would report nothing for ever
ns <- asNamespace("frmtmb")
reg <- function(nm) {
  e <- get(nm, envir = ns, inherits = FALSE)
  stopifnot(is.environment(e))
  e
}
snap <- function() {
  list(
    `compat features` = reg("frmtmb_compat_contrib")$features,
    `compat rules` = reg("frmtmb_compat_contrib")$rules,
    `aterm registry` = reg("frmtmb_aterm_registry")$reg,
    `covstruct registry` = reg("frmtmb_covstruct_registry")$reg,
    # the residual correlation registry, which a frmtmb older than its
    # seam does not have (a "before" arm on the release library)
    `autocor registry` = if (exists("frmtmb_autocor_registry", envir = ns,
                                    inherits = FALSE)) {
      reg("frmtmb_autocor_registry")$reg
    },
    `frame checks` = reg("frmtmb_frame_checks")$fns,
    `prior defaults` = reg("frmtmb_prior_defaults")$providers,
    # a fit gives a notice legitimately, so only a REMOVED tag is a leak
    notices = ls(reg("frmtmb_notice_state"), all.names = TRUE),
    options = options(),
    envvars = as.list(Sys.getenv()),
    `search path` = search(),
    globals = setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed"),
    namespaces = loadedNamespaces())
}
leaks <- function(f, before, after) {
  # A namespace cannot be unloaded safely, so a new one is information,
  # not a leak; it is printed so that a registration it made can be
  # traced to the file that loaded it. An option a namespace adds when
  # it loads is the same: only a CHANGED or removed option is a leak.
  ns_new <- setdiff(after$namespaces, before$namespaces)
  if (length(ns_new)) {
    cat("LOAD ", f, " ", paste(sort(ns_new), collapse = " "), "\n", sep = "")
  }
  added <- setdiff(names(after$options), names(before$options))
  if (length(added)) {
    cat("OPTADD ", f, " ", paste(added, collapse = " "), "\n", sep = "")
  }
  after$options <- after$options[names(after$options) %in%
                                   names(before$options)]
  gone <- setdiff(before$notices, after$notices)
  if (length(gone)) {
    cat("LEAK ", f, " notices: removed ", paste(gone, collapse = " "),
        "\n", sep = "")
  }
  for (k in setdiff(names(before), c("namespaces", "notices"))) {
    b <- before[[k]]
    x <- after[[k]]
    if (identical(b, x)) next
    what <- if (is.list(b) && !is.null(names(b)) && k %in%
                  c("options", "envvars", "aterm registry")) {
      nm <- union(names(b), names(x))
      nm[!vapply(nm, function(n) identical(b[[n]], x[[n]]), TRUE)]
    } else if (is.character(b)) {
      c(paste0("+", setdiff(x, b)), paste0("-", setdiff(b, x)))
    } else {
      paste0(length(b), " -> ", length(x), " entries")
    }
    # an extension registers from its .onLoad(), and a package often
    # sets an option or an environment variable as it loads, so a
    # difference in a file that loaded a namespace names that load; a
    # LEAK line without "with load of" is the test's own. A load creates
    # no global and attaches nothing, so those two are never tagged
    by <- if (k %in% c("globals", "search path")) character(0) else
      if (k %in% c("options", "envvars")) ns_new else
        grep("^frmtmb", ns_new, value = TRUE)
    # rstan defines a Stan model's Rcpp module classes, and the S4
    # class cache entries they bring, in the global environment when it
    # instantiates the model, whoever calls it: rstan::sampling() of a
    # cached model inside a function, with nothing else loaded, leaves
    # the same names (dev/rel073-rcpp-globals.R). They are rstan's, so
    # they are tagged as a load, on a line of their own; every other
    # global stays the test's
    if (identical(k, "globals")) {
      st <- grepl(paste0("^[+](\\.__C__Rcpp_stan_fit4model|",
                         "\\.__C__C[+][+]Object$|\\.__global__$|",
                         "\\.requireCachedGenerics$)"), what)
      if (any(st)) {
        cat("LEAK ", f, " ", k, ": ", paste(what[st], collapse = " "),
            " [with load of a Stan model (rstan)]\n", sep = "")
        what <- what[!st]
        if (!length(what) || all(what == "-")) next
      }
    }
    cat("LEAK ", f, " ", k, ": ", paste(what, collapse = " "),
        if (length(by)) {
          paste0(" [with load of ", paste(sort(by), collapse = " "), "]")
        }, "\n", sep = "")
  }
}

state <- new.env()
Lines <- R6::R6Class("Lines", inherit = testthat::Reporter,
  public = list(
    start_file = function(filename) {
      state$file <- basename(filename)
      state$before <- snap()
      state$t0 <- proc.time()[["elapsed"]]
      cat("FILE ", state$file, "\n", sep = "")
    },
    end_file = function() {
      leaks(state$file, state$before, snap())
      cat("TIME ", state$file, " ",
          round(proc.time()[["elapsed"]] - state$t0, 1), "s\n", sep = "")
    }))

res <- testthat:::test_files(
  test_dir = dir, test_package = p, test_paths = paths,
  load_helpers = TRUE, reporter = Lines$new(), stop_on_failure = FALSE,
  load_package = "installed")

d <- as.data.frame(res)
for (t in res) {
  for (x in t$results) {
    k <- if (inherits(x, "expectation_failure")) "FAIL" else
      if (inherits(x, "expectation_error")) "ERROR" else
        if (inherits(x, "expectation_skip")) "SKIP" else
          if (inherits(x, "expectation_warning")) "WARN" else NA
    if (!is.na(k)) {
      cat("DETAIL ", k, " ", basename(t$file), " [", t$test, "] ",
          gsub("\n", " | ", conditionMessage(x)), "\n", sep = "")
    }
  }
}
# per file in the run's order, in run-tests.R's format
for (f in unique(d$file)) {
  r <- d[d$file == f, , drop = FALSE]
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
# a file that sourced to nothing has no row, and a count that misses it
# is a count of the files that ran, not of the suite
ran <- length(unique(basename(d$file)))
cat("SESSION ", p, " ", ord, " files=", ran, "/", length(paths),
    " pass=", sum(d$passed), " fail=", sum(d$failed), " err=",
    sum(d$error), " skip=", sum(d$skipped), " warn=", sum(d$warning),
    "\n", sep = "")
})
