# Rewrite the no-simulator test of core's test-simulate-density.R for a
# registry without cox(), and drop the cox() test, which moved to
# frmtmb.survival. Run once.
f <- "tests/testthat/test-simulate-density.R"
x <- readLines(f)
a <- grep("^test_that[(]\"a family with no simulator is refused by name at both entry points\",", x)
b <- grep("^test_that[(]\"the deferred categorical[(][)] resolves to a family that simulates\"", x)
stopifnot(length(a) == 1L, length(b) == 1L)
new <- c(
"test_that(\"a family with no simulator is refused by name at both entry points\",",
"          {",
"  # every registry family has a simulator: cox(), the one deliberate",
"  # refusal, is in frmtmb.survival, whose test-cox.R checks it at both",
"  # entry points and checks that it says why",
"  expect_identical(sim_gap_families(), character(0))",
"  # the guard on a family of this file's own: an exponential with mean",
"  # mu and no simulator",
"  tnosim <- frmtmb_family(",
"    \"tnosim\", dpars = \"mu\", links = list(mu = \"log\"),",
"    lpdf = function(y, dpars, aterms) {",
"      -log(dpars[[\"mu\"]]) - y / dpars[[\"mu\"]]",
"    })",
"  set.seed(11)",
"  n <- 60L",
"  dd <- data.frame(x = rep(c(0, 1), length.out = n))",
"  dd$y <- stats::rexp(n, rate = exp(-0.2 + 0.5 * dd$x))",
"  nm <- \"tnosim\"",
"  expect_error(",
"    frmtmb::frm_simulate(y ~ x, dd, family = tnosim,",
"                         newparams = list(b_Intercept = 0, b_x = 0.5),",
"                         nsim = 1L),",
"    paste0(\"family '\", nm, \"' has no simulator\"), fixed = TRUE)",
"  fit <- frmtmb::frm(y ~ x, family = tnosim, data = dd)",
"  expect_error(stats::simulate(fit, nsim = 1L),",
"               paste0(\"family '\", nm, \"' has no simulator\"), fixed = TRUE)",
"})",
"")
x <- c(x[1:(a - 1L)], new, x[b:length(x)])
writeLines(x, f)
