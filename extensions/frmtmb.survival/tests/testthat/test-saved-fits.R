# A saved cox() or royston_parmar() fit read back where this package is
# not installed (review finding 2 of frmtmb's
# dev/reviews/2026-10-23-splitsurv.md). readRDS() then puts the global
# environment in place of this package's namespace, and every internal
# the family's closures call is "not found". frmtmb must stop first, by
# this package's name. The child process sees a library holding a copy
# of the installed frmtmb and none of the libraries that hold this
# package.

saved_child <- function(lib, code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste0(deparse(lib),
                                                collapse = "")),
               code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}

test_that("a saved fit without this package names it, in every method", {
  skip_on_cran()
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  skip_if(is.null(e$bc), "flexsurv's bc data")
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  fr <- frm(bf(recyrs | cens(censored) ~ group),
            family = royston_parmar(df = 2), data = bc)
  set.seed(1)
  dd <- data.frame(x = stats::rnorm(150))
  dd$time <- stats::rexp(150, exp(-0.5 + 0.7 * dd$x))
  dd$out <- stats::rbinom(150, 1, 0.3)
  fc <- frm(bf(time | cens(out) ~ x), family = cox(df = 6), data = dd)
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(list(fr = fr, fc = fc, bc = bc), rds)
  # a library with frmtmb alone, ahead of every library that does not
  # hold this package
  lib <- tempfile("lib")
  dir.create(lib)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  file.copy(system.file(package = "frmtmb"), lib, recursive = TRUE)
  rest <- Filter(function(l) {
    !dir.exists(file.path(l, "frmtmb.survival"))
  }, .libPaths())
  out <- saved_child(c(lib, rest), c(
    "suppressMessages(library(frmtmb))",
    "cat('SURVIVAL ', nzchar(system.file(package = 'frmtmb.survival')),",
    "    '\\n')",
    sprintf("s <- suppressWarnings(readRDS(%s))", deparse(rds)),
    "g <- data.frame(recyrs = c(1, 3),",
    "                group = factor('Good', levels = levels(s$bc$group)))",
    "run <- function(lab, expr) {",
    "  m <- tryCatch({ force(expr); 'NO ERROR' },",
    "                error = function(e) conditionMessage(e))",
    "  cat('CALL ', lab, ' :: ', gsub('\\n', ' ', m), '\\n', sep = '')",
    "}",
    "for (nm in c('fc', 'fr')) {",
    "  f <- s[[nm]]",
    "  run(paste(nm, 'summary'), summary(f))",
    "  run(paste(nm, 'vcov'), vcov(f))",
    "  run(paste(nm, 'fitted'), fitted(f))",
    "  run(paste(nm, 'simulate'), simulate(f, nsim = 1))",
    "  run(paste(nm, 'update'), update(f, . ~ 1))",
    "}",
    "run('fr frm_curve', frm_curve(s$fr, newdata = g, dpar = 'mu',",
    "                              simultaneous = FALSE))",
    "cat('LOGLIK ', sprintf('%.17g', as.numeric(logLik(s$fc))), '\\n')",
    "cat('CHILDOK\\n')"))
  expect_true(any(out == "CHILDOK"), info = paste(out, collapse = "\n"))
  skip_if(any(grepl("^SURVIVAL +TRUE", out)),
          "this package is visible to the child; nothing to check")
  calls <- grep("^CALL ", out, value = TRUE)
  expect_length(calls, 11L)
  # cox() states its own reason for having no simulator before any of
  # its functions is called, so that call needs nothing of this package
  sim <- grepl("^CALL fc simulate", calls)
  expect_identical(sum(sim), 1L)
  expect_match(calls[sim], "family 'cox' has no simulator yet",
               fixed = TRUE)
  calls <- calls[!sim]
  for (cl in calls) {
    expect_match(cl, paste0("is from the extension package ",
                            "frmtmb.survival, which is not installed"),
                 fixed = TRUE, label = cl)
    expect_false(grepl("could not find function", cl, fixed = TRUE),
                 label = cl)
  }
  expect_match(calls[1L], "The cox family of this model", fixed = TRUE)
  expect_match(calls[length(calls)],
               "The royston_parmar family of this model", fixed = TRUE)
  # the guard absent: what needs no function of the family still works,
  # and gives the number this session gives
  ll <- as.numeric(sub("^LOGLIK +", "", grep("^LOGLIK ", out,
                                             value = TRUE)))
  expect_equal(ll, as.numeric(logLik(fc)), tolerance = 1e-12)
})
