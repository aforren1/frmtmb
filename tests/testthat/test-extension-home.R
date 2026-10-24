# A model part from an extension package that this session cannot load
# (review finding 2 of dev/reviews/2026-10-23-splitsurv.md). A saved fit
# keeps its family's closures, and readRDS() without the package that
# built them leaves every internal they call "not found". frmtmb stops
# first, by the package's name, through require_model_packages(). Here
# the family is one this file builds and then says a test package owns,
# a package no session has; frmtmb.survival's test-saved-fits.R does the
# same with real cox() and royston_parmar() fits in a child process.

absent_pkg <- "frmtmbabsenttestpkg"

owned_by <- function(fam, pkg) {
  attr(fam, "frmtmb_package") <- pkg
  fam
}

eh_data <- function() {
  set.seed(4)
  d <- data.frame(x = stats::runif(80))
  d$y <- 1 + 2 * d$x + stats::rnorm(80, 0, 0.3)
  d
}

test_that("require_extension() passes core and a loadable package", {
  expect_true(require_extension("frmtmb", "x"))
  expect_true(require_extension("stats", "x"))
  expect_true(require_extension(NULL, "x"))
  e <- tryCatch(require_extension(absent_pkg, "The tfam family of this model"),
                error = identity)
  expect_s3_class(e, "frmtmb_extension_missing")
  expect_s3_class(e, "frmtmb_error")
  expect_match(conditionMessage(e), paste0(
    "The tfam family of this model is from the extension package ",
    absent_pkg, ", which is not installed"), fixed = TRUE)
  expect_match(conditionMessage(e), paste0("library(", absent_pkg, ")"),
               fixed = TRUE)
})

test_that("a family an older frmtmb built, now an extension's, is refused", {
  # A core family under a moved family's name stands for a cox() fit
  # saved by frmtmb 0.73.0: frmtmb is its package and its density is in
  # frmtmb's namespace. frmtmb's dev/rel074-oldcox.R reads a real one
  fam <- owned_by(as_frmtmb_family(gaussian()), "frmtmb")
  old <- fam
  old[["family"]] <- "cox"
  e <- tryCatch(require_model_packages(old), error = identity)
  expect_s3_class(e, "frmtmb_extension_missing")
  expect_match(conditionMessage(e),
               "The cox family of this model was built by an older frmtmb",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               "Refit the model with library(frmtmb.survival) loaded",
               fixed = TRUE)
  # the guard absent: the core family itself; the same name recorded as
  # the extension's (a frmtmb.eam wiener() fit records frmtmb.eam); a
  # user's family of that name, whose density is not frmtmb's
  expect_null(old_core_family_home(fam))
  expect_true(require_model_packages(fam))
  ext <- old
  ext[["family"]] <- "wiener"
  expect_identical(old_core_family_home(ext), "frmtmb.eam")
  expect_null(old_core_family_home(owned_by(ext, "frmtmb.eam")))
  usr <- old
  lp <- usr[["lpdf"]]
  environment(lp) <- globalenv()
  usr[["lpdf"]] <- lp
  expect_null(old_core_family_home(usr))
})

test_that("a family from an absent package is refused at the fit", {
  d <- eh_data()
  fam <- owned_by(as_frmtmb_family(gaussian()), absent_pkg)
  expect_identical(frm_family_package(fam), absent_pkg)
  expect_error(frm(y ~ x, family = fam, data = d),
               paste0("gaussian family of this model is from the ",
                      "extension package ", absent_pkg), fixed = TRUE)
  # a mixture records its components' packages
  mx <- mixture(gaussian(), fam)
  expect_identical(attr(mx, "frmtmb_component_packages"),
                   c("frmtmb", absent_pkg))
  expect_error(frm(y ~ x, family = mx, data = d), absent_pkg, fixed = TRUE)
  # the guard absent: the same family, owned by core, fits
  expect_s3_class(frm(y ~ x, family = owned_by(fam, "frmtmb"), data = d),
                  "frmtmb_fit")
})

test_that("each method that needs the family's functions names the package", {
  d <- eh_data()
  fit <- frm(y ~ s(x, k = 5), family = gaussian(), data = d)
  bad <- fit
  bad$spec$responses[[1L]]$family <-
    owned_by(bad$spec$responses[[1L]]$family, absent_pkg)
  bad$frame$spec <- bad$spec
  # a fresh cache, as readRDS() gives
  bad$cache <- new.env()
  g <- data.frame(x = seq(0.1, 0.9, length.out = 5))
  # fitted() of a gaussian reads mu through core's own link and so needs
  # nothing of the family's, and update() evaluates gaussian() afresh
  # from the call; frmtmb.survival's test-saved-fits.R checks both on
  # families whose functions are the extension's
  calls <- list(
    summary = function() summary(bad),
    vcov = function() vcov(bad),
    simulate = function() simulate(bad, nsim = 1),
    frm_curve = function() frm_curve(bad, newdata = g, simultaneous = FALSE))
  for (nm in names(calls)) {
    e <- tryCatch({ calls[[nm]](); NULL }, error = identity)
    expect_s3_class(e, "frmtmb_extension_missing")
    expect_match(conditionMessage(e %||% simpleError("")), absent_pkg,
                 fixed = TRUE, label = nm)
  }
  # the guard absent: the fit as it was answers every one
  expect_s3_class(summary(fit), "summary.frmtmb_fit")
  expect_true(is.matrix(vcov(fit)))
  expect_s3_class(frm_curve(fit, newdata = g, simultaneous = FALSE),
                  "frmtmb_curve")
})
