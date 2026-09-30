# Punch round 1: frmtmb.sample test (S01, allow_warnings). Record.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/extensions/frmtmb.sample/tests/testthat/test-formrobust-draws.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())",
"skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

# the sampler's diagnostics on these short chains, and nothing else
fd_ess <- c(\"Effective Samples Size\", \"R-hat\", \"Rhat\", \"lp__\")")
rep1("      ds <- suppressWarnings(suppressMessages(
        frm_sample(bf(y ~ x + ar(t, g, p = 1)), family = gaussian(),
                   data = d, chains = 1, iter = 600, refresh = 0,
                   seed = 3)))",
"      ds <- NULL
      allow_warnings(ds <- suppressMessages(
        frm_sample(bf(y ~ x + ar(t, g, p = 1)), family = gaussian(),
                   data = d, chains = 1, iter = 600, refresh = 0,
                   seed = 3)), fd_ess)")
rep1("        suppressWarnings(suppressMessages(
          frm_sample(bf(f), family = bernoulli(), data = d, chains = 1,
                     iter = 300, refresh = 0, seed = 3)))",
"        out <- NULL
        allow_warnings(out <- suppressMessages(
          frm_sample(bf(f), family = bernoulli(), data = d, chains = 1,
                     iter = 300, refresh = 0, seed = 3)), fd_ess)
        out")
rep1("  i4 <- which(nd$g == \"1\" & nd$t == 4)
  expect_identical(ep[, seq_len(i4)], posterior_epred(cs$ds,
                                                      newdata = full)[
                                                        , seq_len(i4)])",
"  i4 <- which(nd$g == \"1\" & nd$t == 4)
  expect_identical(ep[, seq_len(i4)], posterior_epred(cs$ds,
                                                      newdata = full)[
                                                        , seq_len(i4)])
  # brms fills a missing response with a draw per posterior draw, so the
  # epred draws of the row after it carry the fill's spread, ar * sigma
  # (about 0.55 here), where a row with an observed past carries only
  # the parameters' (about 0.08): brms gives 0.58 against the
  # expected-value fill's 0.10 on these draws
  # (dev/formrobust-rev-log/sample-arma.txt), so the ratio is about 7
  # with the fill and about 1.2 without it
  sd5 <- stats::sd(ep[, i4 + 1L])
  sd6 <- stats::sd(ep[, i4 + 2L])
  expect_gt(sd6 / sd5, 3)")
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\n")
close(con)
