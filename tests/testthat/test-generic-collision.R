# frmtmb must not damage the packages whose generic names it borrows.
#
# Why this file launches CHILD R PROCESSES instead of testing in place.
# The defect is a property of the search path, which is built once per
# session: `library(brms); library(frmtmb)` puts frmtmb's generic above
# brms's, and `UseMethod` then reads frmtmb's method table, which has no
# `brmsfit` entry. The reverse order is not broken at all. A test that
# runs inside the suite's own process sees one order and cannot see the
# other, so each order gets its own process here.
#
# Measured against the 0.55.2 code this replaces (dev/generics-out2,
# dev/generics-summary.R): `library(brms); library(frmtmb)` cost brms
# ALL 27 of the generics the two share, two of them SILENTLY by
# falling into frmtmb's own `.default`, and all 27 of
# `getS3method()`; `library(lme4); library(frmtmb)` cost lme4 all 5 of
# the names they share.
#
# The probe below asks whether the CLASS method is reachable, not
# whether the call errors. An earlier version matched the string "no
# applicable method", and that scored the two silent cases as healthy
# and scored a REGRESSION on `frmtmb_fit` as an improvement, because
# the replacement error did not contain the string it tested for.

# The generics frmtmb exports that brms also has a method for. Kept
# here rather than read out of frmtmb, so that dropping a name from the
# package cannot silently drop it from the test.
brms_shared <- c(
  "as_draws", "as_draws_array", "as_draws_df", "as_draws_list",
  "as_draws_matrix", "as_draws_rvars", "nchains", "ndraws",
  "niterations", "nvariables", "variables", "loo", "loo_compare",
  "waic", "bayes_R2", "prior_summary", "pp_check",
  "conditional_effects", "expose_functions", "hypothesis",
  "posterior_summary", "LOO", "WAIC", "ngrps", "fixef", "ranef",
  "VarCorr"
)

lme4_shared <- c("fixef", "ranef", "VarCorr", "ngrps", "refit")

# Run `code` in a fresh R process on this process's library paths and
# return its stdout. The library paths are passed through because the
# package under test is normally installed somewhere a bare child would
# not look.
run_child <- function(code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  head <- sprintf(".libPaths(%s)",
                  paste0(deparse(.libPaths()), collapse = ""))
  writeLines(c(head, code), f)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(f)), stdout = TRUE, stderr = TRUE))
  paste(out, collapse = "\n")
}

# The probe the children run. `gens` and `cls` are substituted in.
probe_src <- function(attach_code, gens, cls) {
  c(attach_code,
    sprintf("gens <- %s", paste0(deparse(gens), collapse = "")),
    sprintf("obj <- structure(list(), class = %s)", deparse(cls)),
    "bad <- character(); nom <- character()",
    "# UseMethod()'s own lookup: the method table of the namespace",
    "# where the generic that was REACHED is defined.",
    "reach <- function(gen, k) {",
    "  f <- tryCatch(get(gen), error = function(e) NULL)",
    "  if (!is.function(f)) return(NA_character_)",
    "  e <- environment(f)",
    "  if (is.null(e)) e <- baseenv()",
    "  if (!exists('.__S3MethodsTable__.', envir = e,",
    "              inherits = FALSE)) return(NA_character_)",
    "  tb <- get('.__S3MethodsTable__.', envir = e, inherits = FALSE)",
    "  nm <- paste0(gen, '.', k)",
    "  if (exists(nm, envir = tb, inherits = FALSE)) 'yes' else NA",
    "}",
    "for (g in gens) {",
    "  m <- tryCatch(getS3method(g, class(obj), optional = TRUE),",
    "                error = function(e) NULL)",
    "  if (is.null(m)) nom <- c(nom, g)",
    "  if (is.na(reach(g, class(obj)))) bad <- c(bad, g)",
    "}",
    "cat('LOSTDISPATCH:', paste(bad, collapse = ','), '\\n')",
    "cat('LOSTMETHOD:', paste(nom, collapse = ','), '\\n')",
    "cat('CHILDOK\\n')")
}

parse_field <- function(out, key) {
  ln <- grep(paste0("^", key, ":"), strsplit(out, "\n")[[1]], value = TRUE)
  if (!length(ln)) return(NA_character_)
  v <- trimws(sub(paste0("^", key, ":"), "", ln[1]))
  if (!nzchar(v)) character() else strsplit(v, ",")[[1]]
}

test_that("attaching frmtmb after brms leaves every brms method reachable", {
  skip_on_cran()
  skip_if_not_installed("brms")
  out <- run_child(probe_src(
    c("suppressMessages(library(brms))",
      "suppressMessages(library(frmtmb))"),
    brms_shared, "brmsfit"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "LOSTDISPATCH"), character())
  expect_equal(parse_field(out, "LOSTMETHOD"), character())
})

test_that("attaching frmtmb before brms leaves every brms method reachable", {
  skip_on_cran()
  skip_if_not_installed("brms")
  out <- run_child(probe_src(
    c("suppressMessages(library(frmtmb))",
      "suppressMessages(library(brms))"),
    brms_shared, "brmsfit"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "LOSTDISPATCH"), character())
  expect_equal(parse_field(out, "LOSTMETHOD"), character())
})

test_that("brms loaded but never attached is reachable too", {
  # The case the search path cannot rescue: frmtmb is the only thing
  # attached, so its own generic is the one a call finds. Without the
  # load hook in .onLoad this is the one order that stays broken.
  skip_on_cran()
  skip_if_not_installed("brms")
  out <- run_child(probe_src(
    c("suppressMessages(library(frmtmb))",
      "suppressMessages(loadNamespace('brms'))"),
    brms_shared, "brmsfit"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "LOSTDISPATCH"), character())
  expect_equal(parse_field(out, "LOSTMETHOD"), character())
})

test_that("attaching frmtmb after lme4 leaves lme4's own methods reachable", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  out <- run_child(probe_src(
    c("suppressMessages(library(lme4))",
      "suppressMessages(library(frmtmb))"),
    lme4_shared, "merMod"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "LOSTDISPATCH"), character())
  expect_equal(parse_field(out, "LOSTMETHOD"), character())
})

test_that("a shared generic resolves to its owner, not to frmtmb", {
  # The mechanism, stated directly: after the owner is there, the
  # object a user reaches must BE the owner's generic. This is what
  # makes one method table rather than two, and it is the assertion
  # that fails first if the load hook is ever dropped.
  skip_on_cran()
  skip_if_not_installed("posterior")
  skip_if_not_installed("loo")
  out <- run_child(c(
    "suppressMessages(library(posterior))",
    "suppressMessages(library(loo))",
    "suppressMessages(library(frmtmb))",
    "own <- function(g) environmentName(topenv(environment(get(g))))",
    "cat('AS_DRAWS:', own('as_draws_df'), '\\n')",
    "cat('NDRAWS:', own('ndraws'), '\\n')",
    "cat('LOO:', own('loo'), '\\n')",
    "cat('FIXEF:', own('fixef'), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "AS_DRAWS"), "posterior")
  expect_equal(parse_field(out, "NDRAWS"), "posterior")
  expect_equal(parse_field(out, "LOO"), "loo")
  expect_equal(parse_field(out, "FIXEF"), "nlme")
})

test_that("frmtmb loads and works with every optional owner unloadable", {
  # The guard's ABSENT case, built by construction rather than by
  # hoping the machine lacks a package: a directory named for each
  # owner is put FIRST on the library path with a NAMESPACE that
  # stops, so loadNamespace() on it fails however it is reached.
  skip_on_cran()
  lib <- file.path(tempdir(), "frmtmb-absent-owners")
  dir.create(lib, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  owners <- c("posterior", "loo", "rstantools", "bayesplot", "brms",
              "lme4")
  for (p in owners) {
    d <- file.path(lib, p)
    dir.create(file.path(d, "R"), showWarnings = FALSE, recursive = TRUE)
    writeLines(c(paste0("Package: ", p), "Version: 99.0",
                 "Title: Deliberately unusable", "Author: t",
                 "Maintainer: t <t@t.tt>", "Description: t.",
                 "License: GPL-2", "Built: R 4.6.1; ; ; windows"),
               file.path(d, "DESCRIPTION"))
    writeLines("export(nothing)", file.path(d, "NAMESPACE"))
    writeLines(sprintf("stop('%s is deliberately unusable')", p),
               file.path(d, "R", p))
  }
  out <- run_child(c(
    sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),
    "for (p in c('posterior','loo','rstantools','bayesplot','brms',",
    "            'lme4')) {",
    "  ok <- tryCatch({ loadNamespace(p); TRUE }, error = function(e) FALSE)",
    "  if (ok) stop('shadow failed to break ', p)",
    "}",
    "cat('SHADOWED: yes\\n')",
    "suppressMessages(library(frmtmb))",
    "set.seed(20260915)",
    "dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))",
    "dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)",
    "fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)",
    "cat('FIXEF:', paste(names(fixef(fit)), collapse = ','), '\\n')",
    "cat('NGRPS:', ngrps(fit), '\\n')",
    "cat('VARCORR:', length(VarCorr(fit)), '\\n')",
    "cat('REFIT:', class(refit(fit, dd$y))[1], '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "SHADOWED: yes", fixed = TRUE)
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "FIXEF"), c("mu", "sigma"))
  expect_equal(parse_field(out, "VARCORR"), "1")
  expect_equal(parse_field(out, "NGRPS"), "6")
  expect_equal(parse_field(out, "REFIT"), "frmtmb_fit")
})

test_that("an owner unloaded and reloaded does not leave a stale generic", {
  # The construction that broke the round-1 mechanism. A binding
  # SWAPPED once is a hard reference to one namespace instance, so
  # unloading loo and letting something else pull a fresh one left
  # frmtmb dispatching into a dead method table with nothing said.
  # An active binding re-resolves on every access; this pins that.
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("brms")
  out <- run_child(c(
    "suppressMessages(library(loo))",
    "suppressMessages(library(frmtmb))",
    "detach('package:loo'); unloadNamespace('loo')",
    "suppressMessages(loadNamespace('brms'))",
    "g <- get('loo')",
    "cat('ENV:', environmentName(topenv(environment(g))), '\\n')",
    "cat('STALE:', !identical(g, loo::loo), '\\n')",
    "cat('BRMSFIT:', !is.null(getS3method('loo', 'brmsfit',",
    "                                    optional = TRUE)), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "ENV"), "loo")
  expect_equal(parse_field(out, "STALE"), "FALSE")
  expect_equal(parse_field(out, "BRMSFIT"), "TRUE")
})

test_that("the repair survives detach and reattach in both packages", {
  skip_on_cran()
  skip_if_not_installed("brms")
  out <- run_child(c(
    "suppressMessages(library(brms)); suppressMessages(library(frmtmb))",
    "detach('package:frmtmb'); suppressMessages(library(frmtmb))",
    "detach('package:brms'); suppressMessages(library(brms))",
    "cat('LOO:', environmentName(topenv(environment(get('loo')))), '\\n')",
    "cat('CE:', environmentName(topenv(environment(",
    "      get('conditional_effects')))), '\\n')",
    "cat('M:', !is.null(getS3method('conditional_effects', 'brmsfit',",
    "                              optional = TRUE)), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "LOO"), "loo")
  expect_equal(parse_field(out, "CE"), "brms")
  expect_equal(parse_field(out, "M"), "TRUE")
})

test_that("the binding is ACTIVE, and stays active through a re-export", {
  # The mechanism, stated where it can fail. Active-ness propagating
  # through importIntoEnv() is what carries the repair into a package
  # that re-exports the name from frmtmb, which a swapped binding
  # could not do. If this stops being an active binding, the
  # re-export residue comes back and only this assertion says so.
  skip_on_cran()
  out <- run_child(c(
    "suppressMessages(library(frmtmb))",
    "ns <- asNamespace('frmtmb')",
    "cat('NS:', bindingIsActive('loo', ns), '\\n')",
    "cat('ATT:', bindingIsActive('loo',",
    "      as.environment('package:frmtmb')), '\\n')",
    "cat('UNLOCKED:', any(grepl('unlockBinding',",
    "      unlist(lapply(ls(ns), function(n) {",
    "        f <- get0(n, envir = ns); if (is.function(f))",
    "          paste(deparse(body(f)), collapse = ' ') else '' })),",
    "      fixed = TRUE)), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "NS"), "TRUE")
  expect_equal(parse_field(out, "ATT"), "TRUE")
  # and no unlockBinding() anywhere in the namespace, which is the
  # R CMD check NOTE this design exists to avoid
  expect_equal(parse_field(out, "UNLOCKED"), "FALSE")
})

test_that("an owner that exports a NON-generic does not take the name", {
  # frm_owner_generic() checks that the owner's export is a GENERIC,
  # not merely a function. Built rather than assumed: a shadow
  # `bayesplot` whose `pp_check` is a plain function must NOT displace
  # frmtmb's, or frmtmb's own method for its own class becomes
  # unreachable with nothing said.
  #
  # bayesplot is the owner frmtmb registers exactly ONE delayed method
  # on, which is what makes it usable here. Shadowing `posterior` the
  # same way does not test this at all: R resolves
  # `S3method(posterior::as_draws, ...)` against the shadow while
  # LOADING frmtmb, so frmtmb fails to load before any of this runs.
  # That failure is pre-existing and identical on the base build
  # (dev/generics-shadowload.R).
  skip_on_cran()
  lib <- file.path(tempdir(), "frmtmb-nongeneric-owner")
  unlink(lib, recursive = TRUE)
  dir.create(file.path(lib, "bayesplot", "R"), recursive = TRUE,
             showWarnings = FALSE)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  d <- file.path(lib, "bayesplot")
  writeLines(c("Package: bayesplot", "Version: 99.0",
               "Title: A deliberately non-generic pp_check",
               "Author: t", "Maintainer: t <t@t.tt>",
               "Description: t.", "License: GPL-2",
               "Built: R 4.6.1; ; ; windows"),
             file.path(d, "DESCRIPTION"))
  writeLines("export(pp_check)", file.path(d, "NAMESPACE"))
  writeLines("pp_check <- function(object, ...) 42",
             file.path(d, "R", "bayesplot"))
  out <- run_child(c(
    sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),
    "ok <- tryCatch({ loadNamespace('bayesplot'); TRUE },",
    "               error = function(e) FALSE)",
    "cat('SHADOWLOADED:', ok, '\\n')",
    "cat('ISPLAIN:', !any(grepl('UseMethod',",
    "      deparse(body(getExportedValue('bayesplot', 'pp_check'))),",
    "      fixed = TRUE)), '\\n')",
    "suppressMessages(library(frmtmb))",
    "cat('PPCHECK:', environmentName(topenv(environment(",
    "      get('pp_check')))), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  # the guard is only a guard if the shadow really is there and really
  # is not a generic
  expect_equal(parse_field(out, "SHADOWLOADED"), "TRUE")
  expect_equal(parse_field(out, "ISPLAIN"), "TRUE")
  expect_equal(parse_field(out, "PPCHECK"), "frmtmb")
})

test_that("no shared generic carries work in its own body", {
  # The structural guard for a defect that cost 8 assertions and was
  # invisible to a one-file-per-process run.
  #
  # frmtmb's exported generic for any of these names is REPLACEABLE:
  # it resolves to the owner's generic whenever the owner is loaded,
  # and the owner's generic is a bare UseMethod(). So anything frmtmb
  # puts in its own generic runs only when no owner is there.
  # `hypothesis()` used to arm the reserved-name shadowing note in
  # the generic, and under R CMD check, which runs the whole suite in
  # ONE process where an earlier file had loaded brms, the note
  # stopped firing and test-naming-collisions.R lost 8 assertions.
  # Work of that kind belongs in the METHOD.
  #
  # Run in a child with NO owner loaded, so every binding hands back
  # its own fallback rather than the owner's generic.
  skip_on_cran()
  out <- run_child(c(
    "suppressMessages(library(frmtmb))",
    "ns <- asNamespace('frmtmb')",
    "tab <- get('frm_generic_owners', envir = ns)",
    "loaded <- intersect(unique(unlist(tab)), loadedNamespaces())",
    "cat('OWNERSLOADED:', length(loaded), '\\n')",
    "bad <- character()",
    "for (g in names(tab)) {",
    "  f <- tryCatch(get(g), error = function(e) NULL)",
    "  if (!is.function(f)) { bad <- c(bad, g); next }",
    "  b <- body(f)",
    "  if (is.call(b) && identical(b[[1L]], as.name('{')))",
    "    b <- if (length(b) == 2L) b[[2L]] else b",
    "  if (!(is.call(b) && identical(b[[1L]], as.name('UseMethod'))))",
    "    bad <- c(bad, g)",
    "}",
    "cat('NAMES:', length(tab), '\\n')",
    "cat('WITHWORK:', paste(bad, collapse = ','), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  # the guard is only a guard if the fallbacks are what was inspected
  expect_equal(parse_field(out, "OWNERSLOADED"), "0")
  expect_gt(as.integer(parse_field(out, "NAMES")), 20L)
  expect_equal(parse_field(out, "WITHWORK"), character())
})

test_that("every borrowed method frmtmb registers has its owner's twin", {
  # A method frmtmb registers in its OWN table on a borrowed name, say
  # `S3method(loo, frmtmb_fit)`, is reachable only while frmtmb's own
  # generic is the one in use. The moment the owner loads, the exported
  # binding is the owner's generic and the method has to be in the
  # OWNER's table too, which is what `S3method(loo::loo, frmtmb_fit)`
  # puts it. A method added with the first directive and not the second
  # passes every test that runs without the owner loaded and goes
  # unreachable the moment it is loaded. Review found the invariant
  # held, 28 rows with 2 deliberate exceptions, and that nothing
  # asserted it (dev/genrev-r2-pairs.R).
  #
  # The two exceptions are frmtmb's own `.default` methods. Registering
  # those on the owner's generic would REPLACE brms's and loo's own
  # defaults for every class, so they stay in frmtmb's table on purpose.
  ns_file <- parseNamespaceFile("frmtmb", dirname(find.package("frmtmb")))
  m <- ns_file$S3methods
  own <- get("frm_generic_owners", envir = asNamespace("frmtmb"))
  missing_twins <- function(m) {
    local <- which(is.na(m[, 4]) & m[, 1] %in% names(own))
    out <- character()
    for (i in local) {
      for (o in own[[m[i, 1]]]) {
        hit <- m[, 1] == m[i, 1] & m[, 2] == m[i, 2] & m[, 4] %in% o
        if (!any(hit)) out <- c(out, paste0(o, "::", m[i, 1], ".", m[i, 2]))
      }
    }
    sort(out)
  }
  deliberate <- c("brms::posterior_summary.default",
                  "loo::loo_compare.default")
  # the guard is only a guard if it has rows to read
  expect_gt(sum(is.na(m[, 4]) & m[, 1] %in% names(own)), 20L)
  expect_equal(missing_twins(m), deliberate)
  # the inverse, built in: drop ONE real twin and the check must name it
  drop <- which(m[, 1] == "loo" & m[, 2] == "frmtmb_fit" & m[, 4] %in% "loo")
  expect_length(drop, 1L)
  expect_equal(missing_twins(m[-drop, , drop = FALSE]),
               sort(c(deliberate, "loo::loo.frmtmb_fit")))
})

test_that("frm_install_generics() takes a table and refuses a bad one", {
  # frmtmb.sample calls this from its own .onLoad with its own table, so
  # the table is an argument. A malformed one must be refused: a loop
  # over a table with no usable names installs nothing and says nothing,
  # which is the defect back again with a green load.
  bad <- list(list(), list("posterior"), list(a = 1), list(a = character()),
              list(a = NA_character_), setNames(list("loo"), ""))
  for (b in bad) {
    expect_error(frm_install_generics("frmtmb", owners = b),
                 "named list of non-empty character vectors")
  }
  # the inverse: a well-formed table is accepted, and a name frmtmb has
  # already bound is left as it is rather than re-captured
  ns <- asNamespace("frmtmb")
  before <- bindingIsActive("loo", ns)
  expect_true(before)
  expect_identical(frm_install_generics("frmtmb", owners = list(loo = "loo")),
                   "loo")
  expect_true(bindingIsActive("loo", ns))
})
