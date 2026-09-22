# frmtmb.sample must not damage the packages whose generic names it
# borrows.
#
# Why CHILD R PROCESSES. The defect is a property of the search path,
# which is built once per session, so one process sees one load order
# and cannot see another. Each order gets its own process here, the
# same construction as frmtmb's tests/testthat/test-generic-collision.R.
#
# Measured against the tree this replaces (dev/samplegen-check.R, one
# process per order): `library(brms); library(frmtmb.sample)` lost brms
# 28 of the 28 generics this package defines, and so did brms merely
# LOADED, and a package importing frmtmb.sample. The reverse order left
# `rhat()` on a draws object falling SILENTLY into posterior's
# `.default`, and gratia attached after this package did the same to
# `posterior_samples()`.
#
# Every probe below reproduces UseMethod()'s own lookup, the method
# table of the namespace where the reached generic was defined, and
# asks whether the method found there is IDENTICAL to the one that
# should run. It never reads an error string: frmtmb's lane found a
# string-matching counter scoring a silent `.default` as healthy.

# Kept here rather than read from the package, so that dropping a name
# from the package cannot silently drop it from the test.
# log_lik is NOT here: frmtmb defines that generic and the fit refusal,
# and this package re-exports it (R/reexports.R), the way it does loo().
own_generics <- c(
  "as.mcmc", "bayes_factor", "bridge_sampler", "kfold",
  "log_posterior", "loo_moment_match", "loo_subsample", "mcmc_plot",
  "neff_ratio", "nsamples", "nuts_params", "parnames", "post_prob",
  "posterior_epred", "posterior_interval", "posterior_linpred",
  "posterior_predict", "posterior_samples", "pp_mixture",
  "predictive_error", "predictive_interval", "psis", "reloo",
  "restructure", "rhat", "stancode", "standata")

# The owner table, restated for the same reason; a block below asserts
# it is the package's. The first owner of a name is brms's choice where
# two owners disagree (dev/samplegen-audit.R).
owner_table <- list(
  as.mcmc = "coda", bayes_factor = "bridgesampling",
  bridge_sampler = "bridgesampling", post_prob = "bridgesampling",
  kfold = "loo", loo_moment_match = "loo", loo_subsample = "loo",
  psis = "loo", nsamples = "rstantools",
  posterior_epred = "rstantools", posterior_interval = "rstantools",
  posterior_linpred = "rstantools", posterior_predict = "rstantools",
  predictive_error = "rstantools", predictive_interval = "rstantools",
  log_posterior = "bayesplot", neff_ratio = "bayesplot",
  nuts_params = "bayesplot", rhat = c("posterior", "bayesplot"),
  mcmc_plot = "brms", parnames = "brms",
  posterior_samples = c("brms", "gratia"), pp_mixture = "brms",
  reloo = "brms", restructure = "brms", stancode = "brms",
  standata = "brms")

all_owners <- c("posterior", "loo", "rstantools", "bayesplot", "coda",
                "bridgesampling", "gratia", "brms")

test_that("the package installs the owner table this file tests", {
  expect_setequal(names(owner_table), own_generics)
  expect_identical(
    get0("sample_generic_owners", envir = asNamespace("frmtmb.sample"),
         inherits = FALSE),
    owner_table)
})

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

parse_field <- function(out, key) {
  ln <- grep(paste0("^", key, ":"), strsplit(out, "\n")[[1]], value = TRUE)
  if (!length(ln)) return(NA_character_)
  v <- trimws(sub(paste0("^", key, ":"), "", ln[1]))
  if (!nzchar(v)) character() else strsplit(v, ",")[[1]]
}

# A package that imports every one of the names and is loaded, never
# attached. Written as a bare directory, which loadNamespace() accepts.
importer_src <- function(gens) {
  c("implib <- file.path(tempdir(), 'samplegen-importer')",
    "d <- file.path(implib, 'samplegenimp')",
    "dir.create(file.path(d, 'R'), recursive = TRUE, showWarnings = FALSE)",
    "writeLines(c('Package: samplegenimp', 'Version: 0.0.1',",
    "  'Title: t', 'Author: t', 'Maintainer: t <t@t.tt>',",
    "  'Description: t.', 'License: GPL-2',",
    "  'Built: R 4.6.1; ; ; windows'), file.path(d, 'DESCRIPTION'))",
    sprintf("writeLines(c('importFrom(frmtmb.sample, %s)',",
            paste(gens, collapse = ", ")),
    "  'export(nothing)'), file.path(d, 'NAMESPACE'))",
    "writeLines('nothing <- function() NULL',",
    "           file.path(d, 'R', 'samplegenimp'))",
    ".libPaths(c(implib, .libPaths()))")
}

# `where` is the environment the generic is looked up from: the global
# environment for a user, a namespace for a package's own code.
probe_src <- function(setup, where = "globalenv()") {
  c(setup,
    sprintf("gens <- %s", paste0(deparse(own_generics), collapse = "")),
    sprintf("where <- %s", where),
    "tbl <- function(f) get0('.__S3MethodsTable__.',",
    "  envir = environment(f), inherits = FALSE)",
    "reached <- function(g, cls) {",
    "  f <- get0(g, envir = where, mode = 'function')",
    "  if (is.null(f) || is.null(tbl(f))) return(NULL)",
    "  get0(paste0(g, '.', cls), envir = tbl(f), inherits = FALSE)",
    "}",
    "# a method that would run for the class, else the .default that",
    "# would answer it SILENTLY, which counts as lost all the same",
    "runs <- function(g, cls) {",
    "  m <- reached(g, cls)",
    "  if (is.null(m)) m <- reached(g, 'default')",
    "  m",
    "}",
    "lostb <- character(); lostd <- character(); nb <- 0L",
    "sns <- asNamespace('frmtmb.sample')",
    "for (g in gens) {",
    "  mine <- get(paste0(g, '.frmtmb_draws'), envir = sns)",
    "  if (!identical(runs(g, 'frmtmb_draws'), mine)) lostd <- c(lostd, g)",
    "  if (!isNamespaceLoaded('brms')) next",
    "  # the control: brms's method, from the table of brms's own export",
    "  ctl <- get0(paste0(g, '.brmsfit'), inherits = FALSE,",
    "              envir = tbl(getExportedValue('brms', g)))",
    "  if (is.null(ctl)) next",
    "  nb <- nb + 1L",
    "  if (!identical(runs(g, 'brmsfit'), ctl)) lostb <- c(lostb, g)",
    "}",
    "cat('BRMSCONTROL:', nb, '\\n')",
    "cat('LOSTBRMS:', paste(lostb, collapse = ','), '\\n')",
    "cat('LOSTDRAWS:', paste(lostd, collapse = ','), '\\n')",
    "cat('CHILDOK\\n')")
}

expect_clean <- function(out, brms = TRUE) {
  expect_match(out, "CHILDOK", fixed = TRUE)
  if (brms) {
    # the probe is only a probe if it had brms methods to look for
    expect_equal(parse_field(out, "BRMSCONTROL"),
                 as.character(length(own_generics)))
    expect_equal(parse_field(out, "LOSTBRMS"), character())
  }
  expect_equal(parse_field(out, "LOSTDRAWS"), character())
}

test_that("attaching frmtmb.sample after brms loses no brms method", {
  skip_on_cran()
  skip_if_not_installed("brms")
  expect_clean(run_child(probe_src(c(
    "suppressMessages(library(brms))",
    "suppressMessages(library(frmtmb.sample))"))))
})

test_that("attaching frmtmb.sample before brms loses no method either way", {
  # The draws half is what this order tests: brms brings posterior's
  # `rhat`, and a draws method missing from posterior's table fell
  # into posterior's `.default` with nothing said.
  skip_on_cran()
  skip_if_not_installed("brms")
  expect_clean(run_child(probe_src(c(
    "suppressMessages(library(frmtmb.sample))",
    "suppressMessages(library(brms))"))))
})

test_that("brms loaded but never attached loses no brms method", {
  # What a user gets the moment any attached package imports brms.
  skip_on_cran()
  skip_if_not_installed("brms")
  expect_clean(run_child(probe_src(c(
    "suppressMessages(library(frmtmb.sample))",
    "suppressMessages(loadNamespace('brms'))"))))
})

test_that("a package importing frmtmb.sample does not lose brms methods", {
  # Probed from INSIDE the importer's namespace, which is what its own
  # code sees; frmtmb.sample is loaded but never attached.
  skip_on_cran()
  skip_if_not_installed("brms")
  out <- run_child(probe_src(c(
    importer_src(own_generics),
    "suppressMessages(library(brms))",
    "suppressMessages(loadNamespace('samplegenimp'))",
    "cat('ATTACHED:', 'package:frmtmb.sample' %in% search(), '\\n')",
    "cat('IMPACTIVE:', bindingIsActive('posterior_epred',",
    "    parent.env(asNamespace('samplegenimp'))), '\\n')"),
    where = "asNamespace('samplegenimp')"))
  expect_equal(parse_field(out, "ATTACHED"), "FALSE")
  expect_equal(parse_field(out, "IMPACTIVE"), "TRUE")
  expect_clean(out)
})

test_that("the repair survives detach and reattach in both packages", {
  skip_on_cran()
  skip_if_not_installed("brms")
  expect_clean(run_child(probe_src(c(
    "suppressMessages(library(brms))",
    "suppressMessages(library(frmtmb.sample))",
    "detach('package:frmtmb.sample')",
    "suppressMessages(library(frmtmb.sample))",
    "detach('package:brms'); suppressMessages(library(brms))"))))
})

test_that("an owner unloaded and reloaded does not leave a stale generic", {
  skip_on_cran()
  skip_if_not_installed("loo")
  skip_if_not_installed("brms")
  out <- run_child(c(
    "suppressMessages(library(loo))",
    "suppressMessages(library(frmtmb.sample))",
    "detach('package:loo'); unloadNamespace('loo')",
    "cat('UNLOADED:', !isNamespaceLoaded('loo'), '\\n')",
    "suppressMessages(loadNamespace('brms'))",
    "g <- get('psis')",
    "cat('ENV:', environmentName(topenv(environment(g))), '\\n')",
    "cat('STALE:', !identical(g, loo::psis), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "UNLOADED"), "TRUE")
  expect_equal(parse_field(out, "ENV"), "loo")
  expect_equal(parse_field(out, "STALE"), "FALSE")
})

test_that("every other owner's methods survive, in both orders", {
  # brms is one owner of eight. A coda user's as.mcmc() on an mcmc.list
  # or a gratia user's posterior_samples() on a gam is lost the same
  # way. A method counts only if the CONTROL reaches it, the generic the
  # user would find with frmtmb.sample off the search path, because
  # posterior's and bayesplot's rival `rhat` hide each other's methods
  # with no help from this package.
  skip_on_cran()
  owners <- setdiff(all_owners, "brms")
  for (p in owners) skip_if_not_installed(p)
  owner_probe <- c(
    sprintf("own <- %s", paste0(deparse(own_generics), collapse = "")),
    "tbl <- function(f) get0('.__S3MethodsTable__.',",
    "  envir = environment(f), inherits = FALSE)",
    "ctrl <- function(g) {",
    "  for (e in setdiff(search(), 'package:frmtmb.sample')) {",
    "    f <- get0(g, envir = as.environment(e), inherits = FALSE)",
    "    if (is.function(f)) return(f)",
    "  }",
    "  NULL",
    "}",
    "n <- 0L; lost <- character()",
    "for (g in own) {",
    "  cf <- ctrl(g); tf <- get(g)",
    "  if (is.null(cf) || is.null(tbl(cf))) next",
    "  ms <- ls(tbl(cf), all.names = TRUE)",
    "  ms <- ms[startsWith(ms, paste0(g, '.')) &",
    "           !startsWith(ms, 'as.mcmc.list.')]",
    "  for (m in setdiff(ms, paste0(g, '.frmtmb_draws'))) {",
    "    n <- n + 1L",
    "    if (!exists(m, envir = tbl(tf), inherits = FALSE))",
    "      lost <- c(lost, m)",
    "  }",
    "}",
    "cat('OWNERMETHODS:', n, '\\n')",
    "cat('OWNERLOST:', paste(lost, collapse = ','), '\\n')")
  att <- sprintf("suppressMessages(library(%s))", owners)
  for (setup in list(c(att, "suppressMessages(library(frmtmb.sample))"),
                     c("suppressMessages(library(frmtmb.sample))", att))) {
    out <- run_child(probe_src(c(setup, owner_probe)))
    expect_clean(out, brms = FALSE)
    # the guard is only a guard if it found owner methods to look for
    expect_gt(as.integer(parse_field(out, "OWNERMETHODS")), 20L)
    expect_equal(parse_field(out, "OWNERLOST"), character())
  }
})

test_that("the bindings are ACTIVE, in this namespace and not frmtmb's", {
  # frm_install_generics() lives in frmtmb and is called from this
  # package's .onLoad; the bindings must land HERE. Active-ness in the
  # attached environment is what carries the repair to a user.
  skip_on_cran()
  out <- run_child(c(
    "suppressMessages(library(frmtmb.sample))",
    "ns <- asNamespace('frmtmb.sample')",
    sprintf("own <- %s", paste0(deparse(own_generics), collapse = "")),
    "act <- vapply(own, bindingIsActive, NA, env = ns)",
    "cat('NOTACTIVE:', paste(own[!act], collapse = ','), '\\n')",
    "cat('ATT:', bindingIsActive('posterior_epred',",
    "      as.environment('package:frmtmb.sample')), '\\n')",
    "cat('INCORE:', any(vapply(own, exists, NA,",
    "      envir = asNamespace('frmtmb'), inherits = FALSE)), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "NOTACTIVE"), character())
  expect_equal(parse_field(out, "ATT"), "TRUE")
  expect_equal(parse_field(out, "INCORE"), "FALSE")
})

# Eight owners made unusable by construction: a directory named for each
# goes FIRST on the library path with R code that stops, so
# loadNamespace() on it fails however it is reached.
shadow_owners_src <- function() {
  c("shlib <- file.path(tempdir(), 'samplegen-absent-owners')",
    sprintf("owners <- %s", paste0(deparse(all_owners), collapse = "")),
    "for (p in owners) {",
    "  d <- file.path(shlib, p)",
    "  dir.create(file.path(d, 'R'), recursive = TRUE, showWarnings = FALSE)",
    "  writeLines(c(paste0('Package: ', p), 'Version: 99.0',",
    "    'Title: Deliberately unusable', 'Author: t',",
    "    'Maintainer: t <t@t.tt>', 'Description: t.', 'License: GPL-2',",
    "    'Built: R 4.6.1; ; ; windows'), file.path(d, 'DESCRIPTION'))",
    "  writeLines('export(nothing)', file.path(d, 'NAMESPACE'))",
    "  writeLines(sprintf('stop(\"%s is deliberately unusable\")', p),",
    "             file.path(d, 'R', p))",
    "}",
    ".libPaths(c(shlib, .libPaths()))",
    "for (p in owners) {",
    "  ok <- tryCatch({ loadNamespace(p); TRUE }, error = function(e) FALSE)",
    "  if (ok) stop('shadow failed to break ', p)",
    "}",
    "cat('SHADOWED: yes\\n')")
}

fallback_formals <- NULL
test_that("with no owner usable, it loads and its generics are bare", {
  # Two guards on one construction.
  #
  # The absent case: none of the eight owners can load, and the package
  # still loads, dispatches draws to its own methods and refuses by name.
  #
  # The body guard: with no owner loaded every binding hands back this
  # package's own generic, and each must be a bare UseMethod(). While
  # its owner is loaded that generic does not run, so work in its body
  # would stop happening with nothing said. frmtmb found exactly that in
  # `hypothesis()`, and only under R CMD check's single process.
  skip_on_cran()
  out <- run_child(c(
    shadow_owners_src(),
    "suppressMessages(library(frmtmb.sample))",
    sprintf("own <- %s", paste0(deparse(own_generics), collapse = "")),
    "cat('OWNERSLOADED:', paste(intersect(owners, loadedNamespaces()),",
    "    collapse = ','), '\\n')",
    "home <- vapply(own, function(g)",
    "  environmentName(topenv(environment(get(g)))), '')",
    "cat('NOTOURS:', paste(own[home != 'frmtmb.sample'], collapse = ','),",
    "    '\\n')",
    # covr wraps each statement as `if (TRUE) { covr:::count(key); s }`;
    # see the same guard in frmtmb's test-generic-collision.R
    "unwrap <- function(b) {",
    "  repeat {",
    "    if (is.call(b) && identical(b[[1L]], as.name('{')) &&",
    "        length(b) == 2L) { b <- b[[2L]]; next }",
    "    if (is.call(b) && identical(b[[1L]], as.name('if')) &&",
    "        length(b) == 3L && isTRUE(b[[2L]])) {",
    "      i <- b[[3L]]",
    "      if (is.call(i) && identical(i[[1L]], as.name('{')) &&",
    "          length(i) == 3L && is.call(i[[2L]]) &&",
    "          identical(i[[2L]][[1L]], quote(covr:::count))) {",
    "        b <- i[[3L]]; next",
    "      }",
    "    }",
    "    break",
    "  }",
    "  b",
    "}",
    "bad <- character()",
    "for (g in own) {",
    "  b <- unwrap(body(get(g)))",
    "  if (!(is.call(b) && identical(b[[1L]], as.name('UseMethod'))))",
    "    bad <- c(bad, g)",
    "}",
    "cat('WITHWORK:', paste(bad, collapse = ','), '\\n')",
    "for (g in own) cat(sprintf('FORMALS_%s: %s\\n', g,",
    "  gsub(',', ';', paste(deparse(formals(get(g))), collapse = ''))))",
    "set.seed(1)",
    "dd <- data.frame(x = rnorm(40))",
    "dd$y <- rnorm(40)",
    "uf <- frm(bf(y ~ x), family = gaussian(), data = dd,",
    "          dry_run = 'objective')",
    "lab <- c(frmtmb::brms_par_labels(uf), 'lp__')",
    "fd <- structure(list(stanfit = NULL, fit = uf,",
    "  draws = matrix(0, 4, length(lab), dimnames = list(NULL, lab))),",
    "  class = 'frmtmb_draws')",
    "msg <- function(e) tryCatch({ e; 'no error' },",
    "  error = function(c) conditionMessage(c))",
    "cat('STANCODE:', grepl('no Stan program', msg(stancode(fd))), '\\n')",
    "cat('PSAMPLES:', grepl('upgrade path', msg(restructure(fd)),",
    "    fixed = TRUE), '\\n')",
    "cat('BRIDGE:', grepl('integral of the likelihood',",
    "    msg(bridge_sampler(fd))), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "SHADOWED: yes", fixed = TRUE)
  expect_match(out, "CHILDOK", fixed = TRUE)
  # the body guard is only a guard if the fallbacks are what it read
  expect_equal(parse_field(out, "OWNERSLOADED"), character())
  expect_equal(parse_field(out, "NOTOURS"), character())
  expect_equal(parse_field(out, "WITHWORK"), character())
  expect_equal(parse_field(out, "STANCODE"), "TRUE")
  expect_equal(parse_field(out, "PSAMPLES"), "TRUE")
  expect_equal(parse_field(out, "BRIDGE"), "TRUE")
  fallback_formals <<- vapply(own_generics, function(g) {
    paste(parse_field(out, paste0("FORMALS_", g)), collapse = ",")
  }, "")
})

test_that("each generic's formals are its first owner's", {
  # A method must carry every formal of the generic that dispatches to
  # it, and whichever generic dispatches depends on what is loaded. So
  # the fallback here must BE the owner's signature. The first owner in
  # the table is brms's choice where two owners disagree: posterior's
  # `rhat(x, ...)` over bayesplot's `object`, brms's
  # `posterior_samples(x, pars = NA, ...)` over gratia's `model`.
  skip_on_cran()
  skip_if(is.null(fallback_formals),
          "the owner-absent block did not produce the fallback formals")
  tab <- owner_table
  firsts <- unique(vapply(tab, `[[`, "", 1L))
  for (p in firsts) skip_if_not_installed(p)
  out <- run_child(c(
    sprintf("tab <- %s", paste0(deparse(tab), collapse = "")),
    "for (g in names(tab)) {",
    "  f <- getExportedValue(tab[[g]][[1L]], g)",
    "  cat(sprintf('FORMALS_%s: %s\\n', g,",
    "    gsub(',', ';', paste(deparse(formals(f)), collapse = ''))))",
    "}",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  owner_formals <- vapply(own_generics, function(g) {
    paste(parse_field(out, paste0("FORMALS_", g)), collapse = ",")
  }, "")
  differ <- own_generics[owner_formals != fallback_formals]
  expect_equal(differ, character())
})

test_that("an owner exporting a NON-generic does not take the name", {
  # gratia, because this package registers exactly one method on it, so
  # a shadow gratia does not stop the package loading. bbmle's real
  # `parnames` is the same shape: a plain function under an owned name.
  skip_on_cran()
  lib <- file.path(tempdir(), "samplegen-nongeneric-owner")
  unlink(lib, recursive = TRUE)
  d <- file.path(lib, "gratia")
  dir.create(file.path(d, "R"), recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  writeLines(c("Package: gratia", "Version: 99.0",
               "Title: A deliberately non-generic posterior_samples",
               "Author: t", "Maintainer: t <t@t.tt>", "Description: t.",
               "License: GPL-2", "Built: R 4.6.1; ; ; windows"),
             file.path(d, "DESCRIPTION"))
  writeLines("export(posterior_samples)", file.path(d, "NAMESPACE"))
  writeLines("posterior_samples <- function(model, ...) 42",
             file.path(d, "R", "gratia"))
  out <- run_child(c(
    sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),
    "# loadNamespace(), because library() wants a built package",
    "suppressMessages(loadNamespace('gratia'))",
    "cat('ISPLAIN:', !any(grepl('UseMethod',",
    "    deparse(body(gratia::posterior_samples)), fixed = TRUE)), '\\n')",
    "suppressMessages(library(frmtmb.sample))",
    "cat('BRMS:', isNamespaceLoaded('brms'), '\\n')",
    "cat('PS:', environmentName(topenv(environment(",
    "    get('posterior_samples')))), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_equal(parse_field(out, "ISPLAIN"), "TRUE")
  expect_equal(parse_field(out, "BRMS"), "FALSE")
  expect_equal(parse_field(out, "PS"), "frmtmb.sample")
})

test_that("every method on an owned name is in each owner's table", {
  # `S3method(rhat, frmtmb_draws)` reaches only this package's generic.
  # While an owner is loaded the exported name IS the owner's generic, so
  # the method must also be in that owner's table, by
  # `S3method(posterior::rhat, frmtmb_draws)`. Before this invariant was
  # asserted, posterior's `rhat` and gratia's `posterior_samples` were
  # both missing, and a draws object fell into their `.default`.
  ns_file <- parseNamespaceFile("frmtmb.sample",
                                dirname(find.package("frmtmb.sample")))
  m <- ns_file$S3methods
  tab <- owner_table
  missing_twins <- function(m) {
    local <- which(is.na(m[, 4]) & m[, 1] %in% names(tab))
    out <- character()
    for (i in local) {
      for (o in tab[[m[i, 1]]]) {
        hit <- m[, 1] == m[i, 1] & m[, 2] == m[i, 2] & m[, 4] %in% o
        if (!any(hit)) out <- c(out, paste0(o, "::", m[i, 1], ".", m[i, 2]))
      }
    }
    sort(out)
  }
  # the guard is only a guard if it has rows to read
  expect_equal(sum(is.na(m[, 4]) & m[, 1] %in% names(tab)),
               length(own_generics))
  expect_equal(missing_twins(m), character())
  # the inverse, built in: drop ONE real twin and the check must name it
  drop <- which(m[, 1] == "rhat" & m[, 4] %in% "posterior")
  expect_length(drop, 1L)
  expect_equal(missing_twins(m[-drop, , drop = FALSE]),
               "posterior::rhat.frmtmb_draws")
})
