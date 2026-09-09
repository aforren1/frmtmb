# lane tmbstan, probe 06: does GARBAGE satisfy the structural draws
# assertions? This is the claim the lane exists to check: that before
# the static refusal shipped, a CI job could be GREEN while the sampler
# it exercised ignored the model.
#
# Construction. A tmbstan built against StanHeaders >= 2.39 samples a
# standard normal in the UNCONSTRAINED space and then applies Stan's
# constraint transforms, so the returned draws matrix is
# transform(N(0, 1)) column by column, with no data, no likelihood and
# no priors in it. That is synthesized here from a correct run: every
# column is replaced by standard normal noise, and the columns whose
# correct draws are strictly positive (the constrained ones) get
# exp(N(0, 1)) so the object stays well formed. Then the structural
# assertions of test-conditional-effects-draws.R are re-run against it.
#
# SEED 909 for the noise. Model and seeds are those of ce_case() in
# extensions/frmtmb.sample/tests/testthat/test-conditional-effects-draws.R.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-tmbstan/dev/stan-cache")
library(testthat)
library(frmtmb)
library(frmtmb.sample)

band_cols <- c("estimate__", "se__", "lower__", "upper__")
grid_of <- function(df) df[setdiff(names(df), band_cols)]

set.seed(9)
dd <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60),
                 g = factor(rep(1:6, 10)))
dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x + 0.3 * dd$z +
                       stats::rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + z + (1 | g)), family = gaussian(), data = dd)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 2, iter = 500, refresh = 0, seed = 1)))

m <- ds$draws
pos <- apply(m, 2L, function(v) all(v > 0))
cat("draws matrix:", nrow(m), "x", ncol(m), "\n")
cat("columns:", paste(colnames(m), collapse = ", "), "\n")
cat("strictly positive columns:",
    paste(colnames(m)[pos], collapse = ", "), "\n")

set.seed(909)
bad <- m
for (j in seq_len(ncol(bad))) {
  z <- stats::rnorm(nrow(bad))
  bad[, j] <- if (pos[j]) exp(z) else z
}
ds_bad <- ds
ds_bad$draws <- bad
if (!is.null(ds_bad$stanfit)) {
  # the stanfit is only read by diagnostics, not by the structural
  # assertions below; left as is so the object stays valid
}
cat("\nposterior mean, correct run vs garbage:\n")
print(round(rbind(correct = colMeans(m), garbage = colMeans(bad)), 3))

# ---- the structural assertions, verbatim from the test file --------
run <- function(label, expr) {
  n_fail <- 0L
  res <- withCallingHandlers(
    tryCatch({ force(expr); "PASS" },
             expectation_failure = function(e) "FAIL",
             error = function(e) paste0("ERROR: ",
                                        substr(conditionMessage(e), 1, 60))),
    expectation_failure = function(e) {
      n_fail <<- n_fail + 1L
      invokeRestart("muffleCondition")
    })
  cat(sprintf("  %-58s %s (%d failed expectations)\n",
              substr(label, 1, 58), if (n_fail == 0L) res else "FAIL",
              n_fail))
  n_fail == 0L && !startsWith(res, "ERROR")
}

cat("\n== the draws-vs-fit structural block, on GARBAGE draws ==\n")
ok <- c()
ok["columns and grid"] <- run(
  "the draws frame has the fit frame's columns and grid", {
    cd <- conditional_effects(ds_bad, effects = "x", resolution = 8)
    cf <- conditional_effects(fit, effects = "x", resolution = 8)
    expect_identical(names(cd), names(cf))
    expect_identical(names(cd$x), names(cf$x))
    expect_equal(grid_of(cd$x), grid_of(cf$x), ignore_attr = TRUE)
    expect_identical(attr(cd$x, "effects"), attr(cf$x, "effects"))
    expect_true(all(c("y", "z", "g", "cond__", "effect1__") %in%
                      names(cd$x)))
  })
ok["default effect list"] <- run(
  "the default effect list and its grids match the fit method", {
    cd <- conditional_effects(ds_bad, resolution = 6)
    cf <- conditional_effects(fit, resolution = 6)
    expect_identical(names(cd), names(cf))
    for (nm in names(cd)) {
      expect_equal(grid_of(cd[[nm]]), grid_of(cf[[nm]]),
                   ignore_attr = TRUE)
    }
  })
ok["bands are finite and ordered"] <- run(
  "the bands exist, are finite and are ordered", {
    cd <- conditional_effects(ds_bad, effects = "x", resolution = 8)
    expect_true(all(is.finite(cd$x$estimate__)))
    expect_true(all(cd$x$lower__ <= cd$x$estimate__))
    expect_true(all(cd$x$estimate__ <= cd$x$upper__))
  })
ok["draws surface"] <- run(
  "the draws accessor surface answers on garbage", {
    expect_true(nrow(as.matrix(ds_bad)) > 0)
    expect_true(ncol(as.matrix(ds_bad)) > 0)
    ps <- posterior_summary(ds_bad)
    expect_true(all(is.finite(ps[, "Estimate"])))
    pe <- posterior_epred(ds_bad)
    expect_equal(ncol(pe), nrow(dd))
    pp <- posterior_predict(ds_bad)
    expect_equal(ncol(pp), nrow(dd))
  })

cat("\n== the same block on the CORRECT draws, as a control ==\n")
ok2 <- run("columns and grid, correct draws", {
  cd <- conditional_effects(ds, effects = "x", resolution = 8)
  cf <- conditional_effects(fit, effects = "x", resolution = 8)
  expect_identical(names(cd), names(cf))
  expect_equal(grid_of(cd$x), grid_of(cf$x), ignore_attr = TRUE)
})

cat("\n== what a CORRECTNESS assertion says about the same garbage ==\n")
# the ungated correctness assertions are what caught the defect
# upstream; this is the contrast that makes the structural result mean
# something
run("the drawn curve tracks the ML curve", {
  cd <- conditional_effects(ds_bad, effects = "x", resolution = 8)
  cf <- conditional_effects(fit, effects = "x", resolution = 8)
  expect_lt(max(abs(cd$x$estimate__ - cf$x$estimate__)),
            0.5 * stats::sd(dd$y))
})
run("the drawn curve tracks the ML curve, CORRECT draws", {
  cd <- conditional_effects(ds, effects = "x", resolution = 8)
  cf <- conditional_effects(fit, effects = "x", resolution = 8)
  expect_lt(max(abs(cd$x$estimate__ - cf$x$estimate__)),
            0.5 * stats::sd(dd$y))
})

cat("\n== verdict ==\n")
cat("structural blocks passing on garbage:", sum(ok), "of", length(ok),
    "\n")
