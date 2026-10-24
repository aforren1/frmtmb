# cox() is brms's, in brms's core, and frmtmb.survival's here. Guardrail
# 1 of dev/core-boundary.md: every brms spelling of it resolves when
# that package is loaded and is refused by name, naming the package,
# when it is not. frmtmb.survival's test-cox-home.R runs both arms in
# child processes, so that neither depends on what this session loaded;
# these are the arms this session can see.

test_that("family = \"cox\" names frmtmb.survival when it is not loaded", {
  skip_if(isNamespaceLoaded("frmtmb.survival"))
  d <- data.frame(y = c(0.5, 0.7, 0.9), x = 1:3)
  expect_error(frm(y ~ x, family = "cox", data = d),
               "extension package frmtmb.survival, which is not loaded")
  expect_error(brmsfamily("cox"), "frmtmb.survival")
  expect_error(brmsfamily("Cox", link = "identity"), "frmtmb.survival")
  # and the core registry holds no cox() of its own to fall back to
  expect_null(frmtmb:::family_registry[["cox"]])
  expect_identical(frmtmb:::extension_family_homes[["cox"]],
                   "frmtmb.survival")
})

test_that("family = \"cox\" is frmtmb.survival's cox() once that is loaded", {
  # the case the refusal above must not fire on. The name is held in a
  # variable so that R CMD check does not read this as a dependency
  sv <- "frmtmb.survival"
  skip_if_not(requireNamespace(sv, quietly = TRUE))
  expect_identical(frmtmb:::family_ctor("cox"), getExportedValue(sv, "cox"))
  expect_identical(brmsfamily("cox")$family, "cox")
  expect_identical(brmsfamily("cox", link = "identity")$link, "identity")
})
