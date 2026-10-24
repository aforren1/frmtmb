# Guardrail 1 of frmtmb's dev/core-boundary.md for cox(): brms has it in
# its core and frmtmb has it here, so every brms spelling must resolve
# with this package loaded and stop by name, naming this package,
# without it. Both arms run in a child process, so that neither depends
# on what this session has loaded; frmtmb's own test-cox-home.R checks
# the arm its session can see.

home_child <- function(code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  head <- sprintf(".libPaths(%s)",
                  paste0(deparse(.libPaths()), collapse = ""))
  writeLines(c(head, code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}

test_that("without this package, every spelling names it", {
  skip_on_cran()
  out <- home_child(c(
    "suppressMessages(library(frmtmb))",
    "d <- data.frame(y = c(0.5, 0.7, 0.9, 1.3), x = 1:4)",
    "msg <- function(expr) tryCatch({ expr; 'NO ERROR' },",
    "  error = function(e) conditionMessage(e))",
    "cat('STRING ', msg(frm(y ~ x, family = 'cox', data = d)), '\\n')",
    "cat('BRMSFAMILY ', msg(brmsfamily('cox')), '\\n')",
    "cat('LOADED ', isNamespaceLoaded('frmtmb.survival'), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_true(any(out == "CHILDOK"), info = paste(out, collapse = "\n"))
  line <- function(tag) grep(paste0("^", tag, " "), out, value = TRUE)
  want <- "The cox family is in the extension package frmtmb.survival"
  expect_match(line("STRING"), want, fixed = TRUE)
  expect_match(line("STRING"), "library(frmtmb.survival)", fixed = TRUE)
  expect_match(line("BRMSFAMILY"), want, fixed = TRUE)
  # the refusal is the absent case: nothing in the child loaded us
  expect_match(line("LOADED"), "FALSE", fixed = TRUE)
})

test_that("with this package, every spelling is this package's cox()", {
  skip_on_cran()
  out <- home_child(c(
    "suppressMessages(library(frmtmb.survival))",
    "set.seed(1)",
    "d <- data.frame(x = rnorm(120))",
    "d$y <- rexp(120, exp(-0.5 + 0.7 * d$x))",
    "a <- frm(y ~ x, family = 'cox', data = d)",
    "b <- frm(y ~ x, family = brmsfamily('cox'), data = d)",
    "c0 <- frm(y ~ x, family = cox(), data = d)",
    "cat('FAMILY ', family(a)$family, family(b)$family,",
    "    family(c0)$family, '\\n')",
    "cat('SAME ', identical(logLik(a), logLik(c0)),",
    "    identical(logLik(b), logLik(c0)), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_true(any(out == "CHILDOK"), info = paste(out, collapse = "\n"))
  expect_match(grep("^FAMILY ", out, value = TRUE), "cox cox cox",
               fixed = TRUE)
  expect_match(grep("^SAME ", out, value = TRUE), "TRUE TRUE",
               fixed = TRUE)
})

test_that("in this session, the string and brmsfamily() reach cox()", {
  expect_identical(brmsfamily("cox")$family, "cox")
  expect_identical(brmsfamily("cox", link = "identity")$link, "identity")
  set.seed(2)
  d <- data.frame(x = stats::rnorm(120))
  d$y <- stats::rexp(120, exp(-0.5 + 0.7 * d$x))
  a <- frm(y ~ x, family = "cox", data = d)
  expect_identical(stats::family(a)$family, "cox")
  expect_identical(logLik(a), logLik(frm(y ~ x, family = cox(), data = d)))
})
