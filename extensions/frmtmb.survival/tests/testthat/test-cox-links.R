# cox()'s link set, the three checks frmtmb's test-brms-families.R runs
# over every core family, run here for the one family that left core.

test_that("cox()'s link set is the one brms has now", {
  skip_if_not_installed("brms")
  live <- get(".family_cox", envir = asNamespace("brms"))()$links
  expect_identical(live, cox_links)
  expect_identical(cox()$link, live[[1L]])
})

test_that("every link brms accepts for cox() constructs", {
  for (lk in cox_links) {
    expect_identical(cox(link = lk)$link, lk, label = lk)
    expect_identical(brmsfamily("cox", link = lk)$link, lk, label = lk)
  }
  # unquoted, as brms and stats::family() take it
  expect_identical(cox(identity)$link, "identity")
})

test_that("link and linkinv agree with links on cox()", {
  # cox()'s four assertions of frmtmb's "link and link_<dpar> agree with
  # links on every registry family", which loops over core's registry
  f <- cox()
  x <- unclass(f)
  expect_true(is.character(f$link) && length(f$link) == 1L &&
                !is.na(f$link))
  # read from links, never stored on the list
  expect_false(any(grepl("^link($|fun$|inv$|_)", names(x))))
  expect_identical(f$link, x$links$mu$name)
  expect_identical(f$linkinv, x$links$mu$linkinv)
})

test_that("every roster link brms refuses for cox() is refused by name", {
  # frmtmb's link roster, names(frmtmb_links), written out because it
  # is not exported; frmtmb's test-brms-families.R runs the same loop
  # over it for every core family
  roster <- c("identity", "log", "logit", "cloglog", "probit",
              "probit_approx", "cauchit", "softit", "inverse", "logm1",
              "log1p", "tan_half", "power12", "softplus", "squareplus",
              "sqrt", "1/mu^2")
  refused <- 0L
  for (lk in setdiff(roster, cox_links)) {
    expect_error(cox(link = lk),
                 paste0("'", lk, "' is not a supported link for family 'cox'"),
                 fixed = TRUE, label = lk)
    refused <- refused + 1L
  }
  expect_identical(refused, length(setdiff(roster, cox_links)))
  expect_gt(refused, 5L)
  expect_error(cox(link = c("log", "identity")), "single link name")
})
