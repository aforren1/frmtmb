# Reproduce the autotest run recorded in dev/autotest-triage.md.
#
#   Rscript dev/autotest-run.R
#
# The script does three things: it patches three Windows bugs in
# typetracer and one robustness bug in autotest (none of them a frmtmb
# defect, all described in dev/autotest-triage.md), it marks the tests
# that frmtmb declines with the reason, and it runs the rest.
#
# The suppression mechanism is autotest's own: autotest_package() with
# test = FALSE returns one row per test it would run, each with a `test`
# flag; setting a flag to FALSE and handing the table back as
# `test_data` runs everything except those rows. There is no comment
# pragma equivalent to srr's @srrstatsNA, so the exclusions live here,
# next to their reasons, rather than in R/.

PKG <- normalizePath(".", winslash = "/")
OUT <- Sys.getenv("AUTOTEST_OUT", file.path(tempdir(), "autotest"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(NOT_CRAN = "true")

# READ THIS BEFORE TRUSTING ANY RESULT.
#
# typetracer installs the source into a temporary library for tracing,
# and then deletes that library. Every mutation autotest performs
# afterwards resolves the package through .libPaths(), which means
# whatever version happens to be installed on the machine. The first
# run of this script against frmtmb silently tested version 0.28.0 from
# the user library while claiming to test the 0.39.0 source, and a
# third of the rows were version skew rather than diagnostics. Nothing
# in the output says which version was used.
#
# The guard is to install THIS source into a private library and put it
# first, so the fallback is correct whatever typetracer deletes. The
# version actually used is printed at the end; check it.
LIB <- Sys.getenv("AUTOTEST_LIB", file.path(tempdir(), "autotest-lib"))
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(file.path(LIB, "frmtmb"))) {
  message("installing the source under test into ", LIB, " ...")
  utils::install.packages(PKG, lib = LIB, repos = NULL, type = "source",
                          INSTALL_opts = c("--no-multiarch",
                                           "--with-keep.source"))
}
.libPaths(c(LIB, .libPaths()))

library(autotest)

# --- upstream patches -------------------------------------------------

patch_upstream <- function() {
  tt <- asNamespace("typetracer")

  # 1. insert_counters_in_tests() pastes the trace and package
  #    directories into R string literals. On Windows those carry
  #    backslashes, so "C:\Users\..." contains the escape "\U" and every
  #    test file of the traced copy fails to parse.
  ins <- function(pkg_dir) {
    test_path <- file.path(pkg_dir, "tests", "testthat")
    if (!dir.exists(test_path)) return(NULL)
    trace_dir <- gsub("\\\\", "/", options("typetracedir")$typetracedir)
    pkg_fs <- gsub("\\\\", "/", pkg_dir)
    for (f in list.files(test_path, pattern = "^test", recursive = TRUE,
                         full.names = TRUE)) {
      p <- parse(f, keep.source = TRUE)
      out <- lapply(seq_along(p), function(i) {
        pp <- parse(text = deparse(p[[i]]), keep.source = TRUE)
        pd <- utils::getParseData(pp, includeText = TRUE)
        st <- which(pd$token == "SYMBOL_FUNCTION_CALL" &
                      pd$text == "test_that")
        if (!length(st)) return(deparse(p[[i]]))
        sc <- which(pd$token == "STR_CONST")
        sc <- sc[which(sc > st)[1]]
        nm <- gsub("\\\"|\\'", "", pd$text[sc])
        br <- which(pd$token == "'{'")
        br <- br[which(br > sc)[1]]
        if (is.na(br)) return(deparse(p[[i]]))
        idx <- pd$line1[br]
        body <- deparse(p[[i]])
        nm <- gsub("\\s+", "_", nm)
        c(body[seq(idx)], "",
          paste0("traces <- list.files (\"", trace_dir,
                 "\", pattern = \"^typetrace_\", full.names = TRUE)"),
          "ntraces <- length (traces)",
          paste0("ftmp <- file.path (\"", pkg_fs, "\", \"tracetest_", nm,
                 ".txt\")"),
          "writeLines (as.character (ntraces), ftmp)", "",
          body[-seq(idx)])
      })
      writeLines(unlist(out), f)
    }
    NULL
  }
  environment(ins) <- tt
  assignInNamespace("insert_counters_in_tests", ins, ns = tt)

  # 2. reload_pkg() uses tempdir() as a REGULAR EXPRESSION, and a
  #    Windows temp path's backslashes read as back references:
  #    "invalid regular expression ... 'Invalid back reference'".
  rl <- function(pkg_name, lib_path) {
    if (is.null(tryCatch(find.package(pkg_name, lib.loc = .libPaths()),
                         error = function(e) NULL))) {
      return(FALSE)
    }
    fpath <- if (grepl(tempdir(), lib_path, fixed = TRUE)) {
      lib_path
    } else {
      tempdir()
    }
    infile <- file.path(fpath, paste0(pkg_name, "-reload.Rout"))
    outfile <- file.path(fpath, paste0(pkg_name, "-reload-out.Rout"))
    cat("library ('", pkg_name, "')\n", file = infile, sep = "")
    res <- system(paste(shQuote(file.path(R.home("bin"), "R")),
                        "CMD BATCH --vanilla --no-timing",
                        shQuote(infile), shQuote(outfile)))
    # upstream deletes the temp library here, which is what sends every
    # post-tracing mutation to whatever is installed in .libPaths().
    # Keeping it means the traced build stays reachable.
    res == 0L
  }
  environment(rl) <- tt
  assignInNamespace("reload_pkg", rl, ns = tt)

  # 3. trace_package_tests() calls testthat::test_package(), which
  #    throws when any test fails. Under the rewriting above some tests
  #    do fail, and one failure would abort the whole run rather than
  #    the test pass alone. autotest mutates only traces whose source is
  #    "examples", so degrading to examples-only costs parameter
  #    coverage and no mutation.
  orig_tests <- get("trace_package_tests", tt)
  assignInNamespace(
    "trace_package_tests",
    function(package, pkg_dir = NULL, pre_installed = FALSE) {
      tryCatch(orig_tests(package, pkg_dir, pre_installed),
               error = function(e) {
                 message("test tracing degraded to examples-only: ",
                         conditionMessage(e))
                 data.frame(test_file = character(0),
                            test_name = character(0),
                            trace_number = integer(0))
               })
    },
    ns = tt
  )

  # 4. autotest's pass_one_rect_as_other() calls
  #      do.call(data.frame, x$params[[x$i]], quote = TRUE)
  #    which needs the parameter it called rectangular to be a list. A
  #    matrix is rectangular and is not a list, so the call dies with
  #    "second argument must be a list" and takes the run with it.
  at <- asNamespace("autotest")
  orig_rect <- get("pass_one_rect_as_other", at)
  assignInNamespace(
    "pass_one_rect_as_other",
    function(x, other = "data.frame", test_data = NULL) {
      tryCatch(orig_rect(x, other, test_data), error = function(e) NULL)
    },
    ns = at
  )
  invisible(NULL)
}

# --- suppressions -----------------------------------------------------
#
# Each entry names a test class, optionally a function and a parameter,
# and the reason frmtmb declines the expectation. The reasons are
# expanded in dev/autotest-triage.md.

suppressions <- list(
  list(
    test_name = "single_char_case",
    why = paste(
      "Link names, family names, dpar names, prior classes and method",
      "names are identifiers, not free text, and R matches identifiers",
      "by case: stats::binomial(link = 'LOGIT') errors for the same",
      "reason. The one case pair the package DOES accept is",
      "confint(method = 'Wald'), kept because brms spells it that way."
    )
  ),
  list(
    test_name = "par_is_documented",
    why = paste(
      "Every parameter reported here is documented in a combined",
      "roxygen block, e.g. @param lb,ub or @param location,scale,df,",
      "which reaches the .Rd as one \\item{lb, ub}{...}. autotest looks",
      "for an \\item naming one parameter and does not split the",
      "combined form. Verified by reading the \\arguments block of every",
      "reported .Rd."
    )
  ),
  list(
    test_name = "subst_int_for_logical",
    why = paste(
      "autotest expects profile = 1L to keep working. frmtmb refuses",
      "it on purpose: these flags used to be read through isTRUE(),",
      "which calls 1L FALSE, so the integer spelling silently selected",
      "the OPPOSITE option. A refusal is the smaller surprise."
    )
  ),
  list(
    test_name = "single_par_as_length_2", fn_name = "set_prior",
    parameter = c("lb", "ub"),
    why = paste(
      "lb and ub are numeric bounds whose 'unset' marker is NA, so",
      "autotest reads their class off the default and calls them",
      "logical. The mutation tests a type they never had."
    )
  ),
  list(
    test_name = "subst_char_for_logical", fn_name = "set_prior",
    parameter = c("lb", "ub"),
    why = "Same misclassification as above: lb and ub are numeric."
  ),
  list(
    test_name = "random_char_string", fn_name = "set_prior",
    parameter = c("coef", "dpar", "group"),
    why = paste(
      "These name a coefficient, a distributional parameter and a",
      "grouping factor of a model that set_prior() has not seen. Any",
      "string is a well-formed name here; the check that it exists",
      "belongs to the fit the prior is applied to, which does refuse an",
      "unmatched one by name."
    )
  ),
  list(
    test_name = "vector_to_list_col", fn_name = "check_custom_family",
    parameter = "y",
    why = paste(
      "y is the response vector a candidate density is evaluated at.",
      "A list column is not a response, and the density is called on it",
      "directly rather than through a model frame."
    )
  ),
  list(
    test_name = "subst_char_for_logical", fn_name = "frm_bootstrap",
    parameter = "re.form",
    why = paste(
      "re.form defaults to NA and is read as logical for that reason.",
      "It takes NA, NULL or a formula, never TRUE or FALSE."
    )
  ),
  list(
    test_name = c("single_par_as_length_2", "single_char_case"),
    fn_name = "custom_family", parameter = c("dpars", "primary_dpars"),
    why = paste(
      "These are character VECTORS by design, one entry per",
      "distributional parameter. The examples happen to pass a single",
      "name, so autotest records them as length-1 arguments."
    )
  ),
  list(
    test_name = c("single_par_as_length_2", "single_char_case",
                  "random_char_string"),
    fn_name = "frmtmb_family", parameter = "sim_refusal",
    why = paste(
      "sim_refusal is the TEXT of a refusal message. Changing its case",
      "or its length changing the output is correct behavior."
    )
  ),
  list(
    test_name = "negate_logical", fn_name = "frm", parameter = "verbose",
    why = paste(
      "Negating verbose = FALSE to TRUE makes the fit narrate its",
      "stages through message(), which is the documented purpose of the",
      "flag. autotest files each message as its own row, so one flip",
      "produces most of the negate_logical rows in the whole run.",
      "verbose is also the one flag deliberately",
      "left unvalidated: verbose_level() maps anything that is not TRUE,",
      "FALSE or an integer level to silent, and test-verbose.R asserts",
      "that."
    )
  ),
  list(
    test_name = "negate_logical", fn_name = "frm_sample",
    parameter = c("REML", "reparameterize"),
    why = paste(
      "The rows are rstan's own sampler diagnostics (low ESS, one",
      "divergent transition, R-hat up to 1.07) from the short chains the",
      "examples run. They are the sampler reporting on a deliberately",
      "small example, not frmtmb conditions."
    )
  )
)

apply_suppressions <- function(x) {
  x$why <- NA_character_
  for (s in suppressions) {
    hit <- x$test_name %in% s$test_name
    if (!is.null(s$fn_name)) hit <- hit & x$fn_name %in% s$fn_name
    if (!is.null(s$parameter)) hit <- hit & x$parameter %in% s$parameter
    hit[is.na(hit)] <- FALSE
    x$test[hit] <- FALSE
    x$why[hit] <- s$why
  }
  x
}

# --- run --------------------------------------------------------------

patch_upstream()

# Two passes, and each one re-traces every example from scratch: the
# listing pass to learn which tests exist, the testing pass to run them.
# On this package that is about forty minutes. The plan is cached so a
# repeat run is one pass; delete autotest-plan.rds after changing an
# example or a signature.
plan_file <- file.path(OUT, "autotest-plan.rds")
if (file.exists(plan_file)) {
  message("reusing the cached plan in ", plan_file)
  flagged <- readRDS(plan_file)
} else {
  message("listing tests (test = FALSE) ...")
  dummy <- autotest_package(PKG, test = FALSE, progress = "none")
  flagged <- apply_suppressions(dummy)
  saveRDS(flagged, plan_file)
}
message("suppressed ", sum(!flagged$test), " of ", nrow(flagged),
        " candidate tests")

message("running tests (test = TRUE) ...")
res <- autotest_package(PKG, test = TRUE, test_data = flagged,
                        progress = "none")
saveRDS(res, file.path(OUT, "autotest-results.rds"))
utils::write.csv(as.data.frame(res),
                 file.path(OUT, "autotest-results.csv"), row.names = FALSE)
message("rows: ", nrow(res))
print(table(res$type))
print(table(res$test_name, res$type))
message("results in ", OUT)

# The line that says whether any of the above is about this source.
message("frmtmb tested: ",
        as.character(utils::packageVersion("frmtmb")),
        " (from ", dirname(find.package("frmtmb")), ")")
