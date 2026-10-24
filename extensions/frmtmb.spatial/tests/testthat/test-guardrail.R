# Guardrail 1 of frmtmb's dev/core-boundary.md for the specials this
# package took out of core: without the package, frm() refuses each by
# name and says how to load it; with it, a correct model fits.
#
# The refusal arm runs in a CHILD process holding frmtmb alone, because
# this session has the package loaded. The fit arm is this session.

spatial_child <- function(code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste0(deparse(.libPaths()),
                                               collapse = "")),
               "suppressMessages(library(frmtmb))", code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}

# The same, with the child's library paths given.
spatial_child_lib <- function(libs, code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste0(deparse(libs), collapse = "")),
               "suppressMessages(library(frmtmb))", code), f)
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", shQuote(f)), stdout = TRUE,
                           stderr = TRUE))
}

test_that("without frmtmb.spatial, car() and spde() are refused by name", {
  skip_on_cran()
  out <- spatial_child(c(
    "W <- matrix(c(0, 1, 1, 0), 2, dimnames = list(c('a', 'b'), c('a', 'b')))",
    "d <- data.frame(y = rnorm(8), x = rnorm(8), g = rep(c('a', 'b'), 4))",
    "for (fo in c('y ~ x + car(W, gr = g)',",
    "             'y ~ x + car(W, gr = g, type = \"icar\") + (1 | g)',",
    "             'y ~ x + (1 | car(W, gr = g))',",
    "             'y ~ x + spde(fem, gr = node)',",
    "             'y ~ x + (1 | spde(fem, gr = node))')) {",
    "  r <- tryCatch({frm(bf(as.formula(fo)), family = gaussian(),",
    "                     data = d); 'FITTED'},",
    "                error = function(e) conditionMessage(e))",
    "  cat('OUT', gsub('\\n', ' ', r), '\\n')",
    "}",
    "cat('LOADED', 'frmtmb.spatial' %in% loadedNamespaces(), '\\n')"))
  res <- sub("^OUT ", "", grep("^OUT ", out, value = TRUE))
  expect_length(res, 5L)
  for (i in seq_along(res)) {
    nm <- if (i <= 3L) "car" else "spde"
    expect_match(res[[i]], paste0("The ", nm, "() special is in the ",
                                  "extension package frmtmb.spatial, which is ",
                                  "not loaded. Install it and run ",
                                  "library(frmtmb.spatial)"), fixed = TRUE)
  }
  # the child really was without the package, so the refusal is not
  # vacuous
  expect_true(any(out == "LOADED FALSE "), info = paste(out, collapse = "\n"))
})

test_that("with frmtmb.spatial, correct car() and spde() models fit", {
  W <- matrix(0, 4, 4, dimnames = list(paste0("L", 1:4), paste0("L", 1:4)))
  W[cbind(1:3, 2:4)] <- 1
  W <- W + t(W)
  set.seed(1)
  d <- data.frame(loc = factor(rep(rownames(W), each = 5)), x = rnorm(20))
  d$y <- 1 + d$x + rnorm(4)[d$loc] + rnorm(20, 0, 0.5)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")),
             family = gaussian(), data = d, data2 = list(W = W))
  expect_s3_class(fit, "frmtmb_fit")
  expect_true(is.finite(as.numeric(logLik(fit))))
  # the alpha = 2 triple of a regular 1-D mesh, as in test-spde.R
  nn <- 12
  C0 <- diag(rep(0.5, nn))
  C0[1, 1] <- C0[nn, nn] <- 0.25
  G <- matrix(0, nn, nn)
  for (i in seq_len(nn - 1)) {
    G[i, i] <- G[i, i] + 2
    G[i + 1, i + 1] <- G[i + 1, i + 1] + 2
    G[i, i + 1] <- G[i + 1, i] <- -2
  }
  fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
  ds <- data.frame(node = rep(seq_len(nn), each = 6), x = rnorm(6 * nn))
  ds$y <- ds$x + sin(ds$node / 2) + rnorm(6 * nn, 0, 0.3)
  fs <- frm(bf(y ~ x + spde(fem, gr = node)), family = gaussian(),
            data = ds, data2 = list(fem = fem))
  expect_true(is.finite(as.numeric(logLik(fs))))
  # and the column named car is still a column
  d$car <- rnorm(20)
  f2 <- frm(bf(y ~ car + x), family = gaussian(), data = d)
  expect_true("car" %in% names(f2$estimates[["beta"]]))
})

test_that("a saved fit read without frmtmb.spatial is refused by name", {
  skip_on_cran()
  # the fits, saved where the package is loaded
  W <- matrix(0, 4, 4, dimnames = list(paste0("L", 1:4), paste0("L", 1:4)))
  W[cbind(1:3, 2:4)] <- 1
  W <- W + t(W)
  set.seed(1)
  d <- data.frame(loc = factor(rep(rownames(W), each = 5)), x = rnorm(20))
  d$y <- 1 + d$x + rnorm(4)[d$loc] + rnorm(20, 0, 0.5)
  fc <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")), family = gaussian(),
            data = d, data2 = list(W = W))
  nn <- 12
  C0 <- diag(rep(0.5, nn))
  C0[1, 1] <- C0[nn, nn] <- 0.25
  G <- matrix(0, nn, nn)
  for (i in seq_len(nn - 1)) {
    G[i, i] <- G[i, i] + 2
    G[i + 1, i + 1] <- G[i + 1, i + 1] + 2
    G[i, i + 1] <- G[i + 1, i] <- -2
  }
  fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
  ds <- data.frame(node = rep(seq_len(nn), each = 6), x = rnorm(6 * nn))
  ds$y <- ds$x + sin(ds$node / 2) + rnorm(6 * nn, 0, 0.3)
  fs <- frm(bf(y ~ x + spde(fem, gr = node)), family = gaussian(),
            data = ds, data2 = list(fem = fem))
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(rds), add = TRUE)
  saveRDS(list(car = fc, spde = fs), rds)
  every_method <- function(libs) {
    out <- spatial_child_lib(libs, c(
      sprintf("fits <- suppressWarnings(readRDS(%s))", deparse(rds)),
      "cat('HAVE', requireNamespace('frmtmb.spatial', quietly = TRUE), '\\n')",
      "calls <- list(print = function(f) capture.output(print(f)),",
      "  summary = summary, vcov = vcov, predict = predict, fitted = fitted,",
      "  confint = confint, VarCorr = VarCorr, ranef = ranef,",
      "  update = function(f) update(f),",
      "  simulate = function(f) simulate(f, nsim = 1),",
      "  confint_varcorr = confint_varcorr)",
      "for (k in names(fits)) for (m in names(calls)) {",
      "  e <- tryCatch({calls[[m]](fits[[k]]); NULL},",
      "                error = function(e) e)",
      "  cat('CALL', k, m, if (is.null(e)) 'OK' else",
      "      paste(inherits(e, 'frmtmb_extension_missing'),",
      "            gsub('\\n', ' ', conditionMessage(e))), '\\n')",
      "}"))
    expect_true(any(out == "HAVE FALSE "), info = paste(out, collapse = "\n"))
    res <- grep("^CALL ", out, value = TRUE)
    expect_length(res, 22L)
    for (r in res) {
      nm <- strsplit(r, " ", fixed = TRUE)[[1L]][2L]
      expect_match(r, paste0("TRUE The ", nm, "() structure"), fixed = TRUE)
      expect_match(r, paste0("is from the extension package frmtmb.spatial, ",
                             "which is not installed or cannot be loaded"),
                   fixed = TRUE)
    }
  }
  # Installed but unloadable: a frmtmb.spatial that is not a valid
  # installed package, in front of every library this session reads.
  broken <- tempfile("spatial-broken-lib")
  dir.create(file.path(broken, "frmtmb.spatial"), recursive = TRUE)
  on.exit(unlink(broken, recursive = TRUE), add = TRUE)
  writeLines(c("Package: frmtmb.spatial", "Version: 0.0.0"),
             file.path(broken, "frmtmb.spatial", "DESCRIPTION"))
  every_method(c(broken, .libPaths()))
  # Not installed: this session's frmtmb, copied, in front of every
  # library without frmtmb.spatial; the copy masks any other frmtmb there.
  # R CMD check puts frmtmb's dependencies in the library that holds
  # frmtmb.spatial, so there the arm cannot be built and is skipped.
  lib <- tempfile("spatial-absent-lib")
  dir.create(lib)
  on.exit(unlink(lib, recursive = TRUE), add = TRUE)
  file.copy(find.package("frmtmb"), lib, recursive = TRUE)
  rest <- Filter(function(l) {
    !dir.exists(file.path(l, "frmtmb.spatial"))
  }, .libPaths())
  probe <- spatial_child_lib(c(lib, rest), "cat('FRMTMB LOADED\\n')")
  skip_if_not(any(probe == "FRMTMB LOADED"),
              paste("frmtmb's dependencies are only in a library that",
                    "also holds frmtmb.spatial"))
  every_method(c(lib, rest))
})

test_that("loading the package registers its structures and loads no mgcv", {
  skip_on_cran()
  # registering checks each name against a fixed list, so it has no
  # reason to load a namespace the user did not ask for
  out <- spatial_child(c(
    "before <- 'mgcv' %in% loadedNamespaces()",
    "suppressMessages(library(frmtmb.spatial))",
    "cat('MGCV', before, 'mgcv' %in% loadedNamespaces(), '\\n')"))
  expect_true(any(out == "MGCV FALSE FALSE "),
              info = paste(out, collapse = "\n"))
})
