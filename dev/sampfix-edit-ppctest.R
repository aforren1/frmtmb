# Lane sampfix nits round: the multi-chain r_eff test for pp_check().
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/tests/testthat/test-ppcheck-loo.R"
s <- gsub("\r", "", readLines(f))
i <- which(s == "test_that(\"the loo_* types refuse where log_lik() refuses, naming both\", {")
stopifnot(length(i) == 1L)
add <- c(
'test_that("the weights use the chain structure only when every draw is used", {',
'  skip_if_not_installed("loo")',
'  dd <- ppc_data()',
'  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,',
'             dry_run = "objective")',
'  ds <- ppc_draws(fit)',
'  # brms\'s r_eff_log_lik(): relative_eff() on the chains when all draws',
'  # are used, on one chain for a subset. Two chains of 100 here; the',
'  # draws are independent, so the two rules differ only through the',
'  # relative efficiency, which is enough to tell them apart',
'  local_mocked_bindings(draws_nchains = function(x) 2L)',
'  local_mocked_bindings(draws_bayesplot_fun = function(nm, what) {',
'    function(y, yrep, lw, ...) lw',
'  })',
'  ll <- log_lik(ds)',
'  rule <- function(ll, cid) {',
'    stats::weights(loo::psis(-ll, r_eff = loo::relative_eff(',
'      exp(ll), chain_id = cid)), log = TRUE)',
'  }',
'  two <- rule(ll, rep(1:2, each = 100L))',
'  one <- rule(ll, rep(1L, 200L))',
'  expect_gt(max(abs(two - one)), 0)',
'  lw <- allow_warnings(suppressMessages(',
'    pp_check(ds, type = "loo_pit_overlay")), ppc_ok_warnings)',
'  expect_identical(lw, two)',
'  ids <- seq(1L, 200L, by = 3L)',
'  lw_sub <- allow_warnings(suppressMessages(',
'    pp_check(ds, type = "loo_pit_overlay", draw_ids = ids)),',
'    ppc_ok_warnings)',
'  expect_identical(lw_sub, rule(ll[ids, , drop = FALSE],',
'                                rep(1L, length(ids))))',
'})',
'')
writeLines(c(s[seq_len(i - 1L)], add, s[i:length(s)]), f)
cat("ok\n")
