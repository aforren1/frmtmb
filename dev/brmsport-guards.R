# The harness's guards, each with the case where the guarded thing is
# ABSENT, so none of them can pass by never firing.
#
#   Rscript dev/brmsport-guards.R > dev/brmsport-log/guards.txt 2>&1
#
# Every block below MUST report exactly the failures named in its
# label; the script stops otherwise. The gate is checked last, in a
# child process with FRMTMB_BRMS_FIT_TESTS unset.
.libPaths(c("C:/Users/adf44/source/r/brmsport-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
Sys.unsetenv("FRMTMB_BRMSPORT_RECORD")
h <- new.env()
sys.source("tests/testthat/helper-brms-suite.R", envir = h)
attach(h, name = "brmsport-helper")


# every case gets a small fit, as the reviewer's constructions did
run <- function(label, want_fail, expr) {
  f <- tempfile(fileext = ".R")
  writeLines(c("test_that('g', {",
               "set.seed(91)",
               "d <- data.frame(y = rnorm(20), x = rnorm(20))",
               "fit <- frmtmb::frm(y ~ x, d)",
               expr, "})"), f)
  r <- as.data.frame(test_file(f, reporter = "silent",
                               env = new.env(parent = h)))
  got <- sum(r$failed) + sum(r$error)
  cat(sprintf("%-62s failures %d (want %d) %s\n", label, got, want_fail,
              if (got == want_fail) "OK" else "WRONG"))
  got == want_fail
}

res <- c(
  run("pass verdict, assertion holds", 0L,
      "brms_port('g:1', 'pass', '', expect_equal(1, 1))"),
  run("pass verdict, assertion fails", 1L,
      "brms_port('g:2', 'pass', '', expect_equal(1, 2))"),
  run("defect verdict, assertion fails", 0L,
      "brms_port('g:3', 'defect', 'r', expect_equal(1, 2))"),
  run("defect verdict, assertion now holds (stale verdict)", 1L,
      "brms_port('g:4', 'defect', 'r', expect_equal(1, 1))"),
  run("pass verdict, expect_error met by a missing function", 1L,
      "brms_port('g:5', 'pass', '', expect_error(no_such_fn_xyz(), 'no'))"),
  run("pass verdict, expect_error met by a real refusal", 0L,
      "brms_port('g:6', 'pass', '', expect_error(stop('nope'), 'no'))"),
  run("pass verdict, bare expect_error met by an argument refusal", 1L,
      paste0("brms_port('g:7', 'pass', '', expect_error(",
             "frmtmb::frm(y ~ 1, data.frame(y = 1:3), bogus_arg = 1)))")),
  run("pass verdict, reads an object whose assignment failed", 1L,
      c("x <- 1",
        "brms_setup('g:8', x <- stop('setup broke'))",
        "brms_port('g:9', 'pass', '', expect_equal(x, 1))")),
  run("pass verdict, object reassigned successfully after a failure", 0L,
      c("brms_setup('g:10', x <- stop('setup broke'))",
        "brms_setup('g:11', x <- 1)",
        "brms_port('g:12', 'pass', '', expect_equal(x, 1))")),
  run("not-run assertion recorded as pass", 1L,
      "brms_port_not_run('g:13', 'pass', '', 'why')"),
  run("hollow verdict, still holds", 0L,
      "brms_port('g:14', 'hollow', 'r', expect_true(is.null(NULL)))"),
  run("hollow verdict, stops holding", 1L,
      "brms_port('g:15', 'hollow', 'r', expect_true(is.null(1)))"),
  # punch round 1: the reviewer's seven hollow constructions, each beside
  # the genuine case its rule must NOT reject
  run("H1 regex expect_error met by an argument-name refusal", 1L,
      "brms_port('h:1', 'pass', '', expect_error(fitted(fit, bogus = 1), 'bogus'))"),
  run("H1 control: the same regex met by a value refusal", 0L,
      "brms_port('h:1c', 'pass', '', expect_error(stop('bogus value'), 'bogus'))"),
  run("H2 bare expect_error met by a 'cannot honor' refusal", 1L,
      "brms_port('h:2', 'pass', '', expect_error(fitted(fit, ndraws = 5)))"),
  run("H2 control: bare expect_error met by a real refusal", 0L,
      "brms_port('h:2c', 'pass', '', expect_error(log(-1:1, base = 'a')))"),
  run("H3 expect_equal on two NULLs through fit$data", 1L,
      "brms_port('h:3', 'pass', '', expect_equal(fit$data$y, fit$data$x))"),
  run("H3 control: expect_equal on a NULL against the literal NULL", 0L,
      "brms_port('h:3c', 'pass', '', expect_equal(list()$x, NULL))"),
  run("H3 control: expect_equal on two equal non-NULL values", 0L,
      "brms_port('h:3d', 'pass', '', expect_equal(nobs(fit), 20))"),
  run("H4 regex met by R's missing-argument error", 1L,
      "brms_port('h:4', 'pass', '', expect_error((function(group) group)(), 'group'))"),
  run("H4 control: the regex met by a designed refusal", 0L,
      "brms_port('h:4c', 'pass', '', expect_error(stop('group is required'), 'group'))"),
  run("H5 stale object through x$a <-", 1L,
      c("x <- list(a = 1)",
        "brms_setup('h:5', x$a <- stop('setup broke'))",
        "brms_port('h:6', 'pass', '', expect_equal(x$a, 1))")),
  run("H5 control: a successful x$a <- leaves x fresh", 0L,
      c("x <- list(a = 1)",
        "brms_setup('h:5c', x$a <- 2)",
        "brms_port('h:6c', 'pass', '', expect_equal(x$a, 2))")),
  run("H5 control: x$a <- after a failed x <- keeps x stale", 1L,
      c("x <- list(a = 1)",
        "brms_setup('h:5d', x <- stop('broke'))",
        "brms_setup('h:5e', x$a <- 2)",
        "brms_port('h:6d', 'pass', '', expect_equal(x$a, 2))")),
  run("H6 a missing function 'of mode function was not found'", 1L,
      "brms_port('h:7', 'pass', '', expect_error(match.fun('no_such_fn_xyz')(), 'no_such'))"),
  # H7 has no rule since punch round 2: is.numeric() on an NA cannot be
  # told from an NA that is the right answer (V3), so row 394 is hollow
  # by hand, and the hand marking is what is guarded
  run("H7/V3 is.numeric() on NA is not rejected by a rule", 0L,
      "brms_port('h:8', 'pass', '', expect_true(is.numeric(NA_real_)))"),
  run("H7 hand-marked hollow: still holds", 0L,
      "brms_port('h:8h', 'hollow', 'r', expect_true(is.numeric(NA_real_)))"),
  run("H7 hand-marked hollow: stops holding", 1L,
      "brms_port('h:8i', 'hollow', 'r', expect_true(is.numeric('a')))"),
  # punch round 2: the recheck's under-reach (U) and over-reach (V) forms
  run("U1 expect_null on a partial $ read", 1L,
      "brms_port('u:1', 'pass', '', expect_null(fit$data$y))"),
  run("U1 control: expect_null on a final absent name", 0L,
      "brms_port('u:1c', 'pass', '', expect_null(list(a = 1)$b))"),
  run("U2 expect_true(is.null()) on a partial $ read", 1L,
      "brms_port('u:2', 'pass', '', expect_true(is.null(fit$data$y)))"),
  run("U5 expect_length on a partial $ read", 1L,
      "brms_port('u:5', 'pass', '', expect_length(fit$data$y, 0))"),
  run("U5 control: expect_length on a real element", 0L,
      "brms_port('u:5c', 'pass', '', expect_length(d$y, 20))"),
  run("U6 argument refusal in other words", 1L,
      "brms_port('u:6', 'pass', '', expect_error(frmtmb::bf(y ~ x, sigma1 = 'sigma2'), 'sigma'))"),
  run("U7 stale through assign()", 1L,
      c("x <- 1", "brms_setup('u:7', assign('x', stop('broke')))",
        "brms_port('u:8', 'pass', '', expect_equal(x, 1))")),
  run("U7 control: a successful assign() is fresh", 0L,
      c("x <- 1", "brms_setup('u:7c', assign('x', 2))",
        "brms_port('u:8c', 'pass', '', expect_equal(x, 2))")),
  run("U8 NULL against a variable holding NULL, partial read", 1L,
      c("e <- NULL",
        "brms_port('u:9', 'pass', '', expect_equal(fit$data$y, e))")),
  run("V1 two genuinely NULL names() pass", 0L,
      "brms_port('v:1', 'pass', '', expect_equal(names(c(1, 2)), names(c(3, 4))))"),
  run("V2 an assertion ABOUT an argument refusal passes", 0L,
      "brms_port('v:2', 'pass', '', expect_error(fitted(fit, bogus = 1), 'has no argument'))"),
  run("O2 own words reading a stale object", 1L,
      c("x <- fit", "brms_setup('o:2s', x <- stop('broke'))",
        "brms_port_own('o:3s', \"Unknown dpar: 'inv'\", '', expect_error(fitted(x, dpar = 'inv'), 'Invalid argument'))")),
  run("O2 control: the same own-words row on a fresh object", 0L,
      "brms_port_own('o:3f', \"Unknown dpar: 'inv'\", '', expect_error(fitted(fit, dpar = 'inv'), 'Invalid argument'))"),
  run("O3 own-words pattern '.' is not specific", 1L,
      "brms_port_own('o:6', '.', '', expect_error(frmtmb::frm(y ~ x, d, family = poisson()), 'brms words'))"),
  run("O3 own-words pattern matching the empty string", 1L,
      "brms_port_own('o:6e', 'x*', '', expect_error(frmtmb::frm(~ x, d), 'brms words'))"),
  run("O3 own-words pattern matching an unrelated refusal only", 1L,
      "brms_port_own('o:6u', 'is not supported', '', expect_error(frmtmb::frm(~ x, d), 'brms words'))"),
  # punch round 2, recorded and not built: a type check on an NA cannot
  # be told from an NA that is the right answer (V3), so neither form is
  # rejected, and a hollow one must be marked by hand as row 394 is
  run("U3 expect_type() on NA is not rejected by a rule", 0L,
      "brms_port('u:3', 'pass', '', expect_type(NA_real_, 'double'))"),
  run("U4 compound type check on NA is not rejected by a rule", 0L,
      c("v <- NA_real_",
        "brms_port('u:4', 'pass', '', expect_true(is.numeric(v) && length(v) == 1))")),
  run("reviewer control: bare expect_error met by 'has no argument'", 1L,
      "brms_port('h:9', 'pass', '', expect_error(vcov(fit, cor = TRUE)))"),
  run("reviewer control: simple stale object", 1L,
      c("x <- 1", "brms_setup('h:10', x <- stop('broke'))",
        "brms_port('h:11', 'pass', '', expect_equal(x, 1))")),
  # the own-words verdict of the user's rule 1
  run("own words: frmtmb's refusal of the same case holds", 0L,
      "brms_port_own('o:1', 'needs a response', '', expect_error(frmtmb::frm(~ x, d), 'Response variable is missing'))"),
  run("own words: a DIFFERENT error does not satisfy it (brm:81)", 1L,
      "brms_port_own('o:2', 'se[(][)]', '', expect_error(frmtmb::frm(y | se(sei) ~ x, d), 'se'))"),
  run("own words: no error at all does not satisfy it", 1L,
      "brms_port_own('o:3', 'needs a response', '', expect_error(frmtmb::frm(y ~ x, d), 'Response variable is missing'))"),
  run("own words: brms's pattern already holds (stale verdict)", 1L,
      "brms_port_own('o:4', 'response', '', expect_error(frmtmb::frm(~ x, d), 'needs a response'))"),
  run("own words: fixed = TRUE in brms's call is dropped", 0L,
      "brms_port_own('o:5', 'needs a response [(]left', '', expect_error(frmtmb::frm(~ x, d), fixed = TRUE, 'Response variable (is) missing'))")
)

# the frmtmb.sample copy of the helper cannot drift from core's: every
# function in one is identical() to the same function in the other,
# checked on the real copy and on a copy with one body changed
helper_same <- function(core_path, copy_path) {
  a <- new.env()
  b <- new.env()
  sys.source(core_path, envir = a, keep.source = FALSE)
  sys.source(copy_path, envir = b, keep.source = FALSE)
  na <- sort(ls(a, all.names = TRUE))
  nb <- sort(ls(b, all.names = TRUE))
  identical(na, nb) && all(vapply(na, function(n) {
    x <- get(n, a)
    y <- get(n, b)
    if (is.function(x)) {
      identical(deparse(x), deparse(y))
    } else if (is.environment(x)) {
      TRUE
    } else {
      identical(x, y)
    }
  }, TRUE))
}
core_h <- "tests/testthat/helper-brms-suite.R"
copy_h <- "extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R"
same_real <- helper_same(core_h, copy_h)
mut <- tempfile(fileext = ".R")
txt <- readLines(copy_h)
i <- grep("^brms_shim_rename_pars <- function", txt)
txt[i] <- "brms_shim_rename_pars <- function(x) NULL"
writeLines(txt, mut)
same_mut <- helper_same(core_h, mut)
cat(sprintf("helper copy identical to core: %s; a mutated copy: %s\n",
            same_real, same_mut))
res <- c(res, same_real, !same_mut)

# The same identity as a TEST (recheck R6): frmtmb.sample's tier file
# test-brms-suite-helper-copy.R must pass on the real tree and fail on a
# tree whose copy has one body changed. Both trees are built in a
# temporary directory with the layout the test's relative path expects.
copy_test <- function(copy_src) {
  root <- tempfile("brmsport-copytree")
  cdir <- file.path(root, "tests", "testthat")
  sdir <- file.path(root, "extensions", "frmtmb.sample", "tests",
                    "testthat")
  dir.create(cdir, recursive = TRUE)
  dir.create(sdir, recursive = TRUE)
  file.copy(core_h, cdir)
  file.copy(copy_src, file.path(sdir, "helper-brms-suite.R"))
  tf <- file.path(sdir, "test-brms-suite-helper-copy.R")
  file.copy(file.path("extensions/frmtmb.sample/tests/testthat",
                      "test-brms-suite-helper-copy.R"), tf)
  old <- Sys.getenv("FRMTMB_BRMS_FIT_TESTS")
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
  on.exit(Sys.setenv(FRMTMB_BRMS_FIT_TESTS = old))
  r <- as.data.frame(test_file(tf, reporter = "silent",
                               env = new.env(parent = h)))
  c(passed = sum(r$passed), failed = sum(r$failed) + sum(r$error),
    skipped = sum(r$skipped))
}
ct_real <- copy_test(copy_h)
ct_mut <- copy_test(mut)
cat(sprintf(paste("helper-copy test on the real copy: pass %d fail %d",
                  "skip %d; on a mutated copy: pass %d fail %d skip %d\n"),
            ct_real[["passed"]], ct_real[["failed"]], ct_real[["skipped"]],
            ct_mut[["passed"]], ct_mut[["failed"]], ct_mut[["skipped"]]))
# the real side must also have asserted something, or a skip reads clean
res <- c(res,
         ct_real[["failed"]] == 0 && ct_real[["skipped"]] == 0 &&
           ct_real[["passed"]] > 0,
         ct_mut[["failed"]] > 0)

# Specificity against REAL refusals (recheck R1): no own-words pattern
# may match the frmtmb message another own-words row caught for its own,
# different case. Read from the record, so run this after
# dev/brmsport-record.sh. The inverse case: a pattern built to be
# unspecific ("response") is seen to match another row's message.
recs <- do.call(rbind, lapply(
  Sys.glob("dev/brmsport-log/rec-*.tsv"), function(f) {
    utils::read.delim(f, header = FALSE, quote = "",
                      colClasses = "character", na.strings = NULL,
                      col.names = c("kind", "pkg", "id", "verdict", "held",
                                    "vacuous", "msg", "caught", "raw"))
  }))
own_rows <- unique(recs[recs$verdict == "own" & recs$held == "TRUE",
                        c("id", "caught")])
own_pat <- utils::read.delim("dev/brmsport-verdicts-own.tsv", quote = "",
                             colClasses = "character", na.strings = NULL)
cross <- function(pats, ids) {
  hits <- character()
  for (k in seq_along(pats)) {
    other <- own_rows$caught[own_rows$id != ids[k]]
    if (any(grepl(pats[k], other))) hits <- c(hits, ids[k])
  }
  hits
}
cross_real <- cross(own_pat$pattern, own_pat$id)
cross_bad <- cross("response", "standata:106")
cat(sprintf(paste("own-words rows with a caught message: %d of %d;",
                  "patterns matching another row's message: %d (%s);",
                  "an unspecific control pattern matches: %d\n"),
            length(unique(own_rows$id)), nrow(own_pat), length(cross_real),
            paste(cross_real, collapse = " "), length(cross_bad)))
res <- c(res, length(unique(own_rows$id)) == nrow(own_pat),
         length(cross_real) == 0L, length(cross_bad) == 1L)


# the gate: with FRMTMB_BRMS_FIT_TESTS unset a generated file skips whole
gate_run <- function(set) {
  gs <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf(".libPaths(c(%s))",
            paste(sprintf("'%s'", .libPaths()), collapse = ", ")),
    "Sys.setenv(NOT_CRAN = 'true')",
    if (set) "Sys.setenv(FRMTMB_BRMS_FIT_TESTS = 'true')" else
      "Sys.unsetenv('FRMTMB_BRMS_FIT_TESTS')",
    "suppressMessages({library(testthat); library(frmtmb)})",
    "r <- as.data.frame(test_file(",
    "  'tests/testthat/test-brms-suite-families.R', package = 'frmtmb',",
    "  env = test_env('frmtmb'), reporter = 'silent'))",
    "writeLines(paste('GATE', sum(r$passed)))"), gs)
  out <- system2(file.path(R.home("bin"), "Rscript.exe"), gs,
                 stdout = TRUE, stderr = TRUE)
  hit <- grep("^GATE [0-9]+$", out, value = TRUE)
  as.integer(substring(hit, 6L))
}
off <- gate_run(FALSE)
on <- gate_run(TRUE)
cat("families file passes with the variable unset:", off, "and set:", on,
    fill = TRUE)
# the unset side alone would also read 0 on a path that ran nothing
gate_ok <- length(off) == 1L && length(on) == 1L && off == 0L && on == 84L
cat(sum(res), "of", length(res), "guards behave;",
    if (gate_ok) "the gate skips" else "THE GATE CHECK FAILED", fill = TRUE)
stopifnot(all(res), gate_ok)
