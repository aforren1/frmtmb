# The test helper local_frmtmb_registries() (helper-registry.R) is what
# keeps a registering test from reaching the files after it, so its own
# failure modes are pinned here: a registry it does not know about, and
# an extension that loads inside its scope.

test_that("the helper knows every environment the namespace holds", {
  # A new registry, or a renamed one, would be saved by nobody and
  # restored by nobody, and a test that registered into it would leak
  # in silence. Every top-level environment is either a registry the
  # helper saves or a cache it leaves alone on purpose.
  ns <- asNamespace("frmtmb")
  envs <- Filter(function(nm) {
    !bindingIsActive(nm, ns) &&
      is.environment(get(nm, envir = ns, inherits = FALSE)) &&
      !nm %in% c(".__NAMESPACE__.", ".__S3MethodsTable__.")
  }, ls(ns, all.names = TRUE))
  expect_setequal(envs, c(frmtmb_test_registries, frmtmb_test_caches))
  # and the helper reads them so that an absent name stops it
  expect_error(frmtmb_test_registry("frmtmb_no_such_registry"))
})

test_that("an extension loaded inside the scope keeps its registrations", {
  skip_on_cran()
  pk <- c("frmtmb.eam", "frmtmb.learn", "frmtmb.sample")
  have <- vapply(pk, function(p) nzchar(system.file(package = p)), TRUE)
  skip_if_not(all(have), paste("not installed:", toString(pk[!have])))
  helper <- normalizePath(testthat::test_path("helper-registry.R"), "/")
  # one child loads the extensions inside the helper's scope, the other
  # outside any scope; the registries must come out the same. learn
  # imports eam, so the scope sees two namespaces appear and has to run
  # eam's registration before learn's. frmtmb.sample's .onLoad() also
  # installs its generics as active bindings, which must survive the
  # second call
  probe <- c(
    "ns <- asNamespace('frmtmb')",
    "st <- function() list(",
    "  features = ns$frmtmb_compat_contrib$features,",
    "  rules = length(ns$frmtmb_compat_contrib$rules),",
    "  aterms = names(ns$frmtmb_aterm_registry$reg),",
    "  checks = length(ns$frmtmb_frame_checks$fns),",
    "  priors = length(ns$frmtmb_prior_defaults$providers))",
    "s <- st()",
    "cat('STATE', gsub('\\n', ' ', paste(deparse(s), collapse = '')),",
    "    '\\n')",
    "r <- tryCatch({",
    "  frm(bf(y | dec(d) ~ x) + gaussian(), data = NULL,",
    "      dry_run = 'spec'); 'parses'",
    "}, error = function(e) 'refused')",
    "cat('DEC', r, '\\n')",
    "sns <- asNamespace('frmtmb.sample')",
    "act <- Filter(function(n) bindingIsActive(n, sns), ls(sns))",
    "ok <- all(vapply(act, function(n) is.function(get(n, sns)), TRUE))",
    "cat('GEN', length(act), ok, '\\n')")
  head <- c(sprintf(".libPaths(%s)",
                    paste0(deparse(.libPaths()), collapse = "")),
            "suppressMessages(library(frmtmb))",
            sprintf("source(%s)", deparse(helper)))
  run <- function(code) {
    f <- tempfile(fileext = ".R")
    on.exit(unlink(f), add = TRUE)
    writeLines(c(head, code, probe), f)
    suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                             c("--vanilla", shQuote(f)), stdout = TRUE,
                             stderr = TRUE))
  }
  scoped <- run(c("f <- function() {",
                  "  local_frmtmb_registries()",
                  "  for (p in c('frmtmb.learn', 'frmtmb.sample')) {",
                  "    loadNamespace(p)",
                  "  }",
                  "}",
                  "f()"))
  plain <- run("for (p in c('frmtmb.learn', 'frmtmb.sample')) loadNamespace(p)")
  st <- function(out) grep("^STATE ", out, value = TRUE)
  expect_length(st(plain), 1L)
  expect_identical(st(scoped), st(plain),
                   info = paste(scoped, collapse = "\n"))
  expect_identical(grep("^DEC ", scoped, value = TRUE), "DEC parses ")
  gen <- function(out) grep("^GEN ", out, value = TRUE)
  expect_identical(gen(scoped), gen(plain))
  expect_match(gen(plain), "^GEN [1-9][0-9]* TRUE")
})
