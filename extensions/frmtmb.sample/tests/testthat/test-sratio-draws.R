# sratio()'s thresholds are an unconstrained vector, as brms declares
# them, so the sampler's `tau_raw` draws ARE the thresholds and a draw
# may have them cross. frmtmb held them as (first threshold, log
# increments) before, which kept every draw ordered.

test_that("sratio threshold draws are the thresholds and may cross", {
  skip_on_cran()
  skip_sampler()
  # the hazard of stopping at category 2 is well below that at 1, so
  # the second threshold lies below the first
  counts <- c(50, 10, 40, 25)
  d <- data.frame(y = rep(seq_along(counts), counts))
  fit <- frm(y ~ 1, data = d, family = sratio())
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 4)))
  fe <- fixef(ds, summary = FALSE)
  raw <- as.matrix(ds, variable = paste0("tau_raw_", 1:3))
  # as.vector(): the draws matrix carries an nchains attribute
  expect_identical(as.vector(fe[, paste0("Intercept[", 1:3, "]")]),
                   as.vector(raw))
  expect_gt(mean(fe[, "Intercept[2]"] < fe[, "Intercept[1]"]), 0.5)
})
