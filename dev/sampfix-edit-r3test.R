# Lane sampfix last round (R3): a later draw that overflows is not a read.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/tests/testthat/test-laplace-draws.R"
s <- gsub("\r", "", readLines(f))
i <- which(s == "test_that(\"a refusal names the function the user called\", {")
stopifnot(length(i) == 1L)
add <- c(
"test_that(\"a later draw that overflows on its own is not taken for a read\", {",
"  # exp(720) is Inf at draw 3 for a reason that has nothing to do with",
"  # the integrated values, and the full draws return it; a read of those",
"  # values comes out NA, never Inf (dev/sampfix-12-fillpaths.R)",
"  set.seed(77)",
"  dd <- data.frame(g = factor(rep(1:8, each = 10)))",
"  dd$x <- stats::rnorm(nrow(dd))",
"  u <- stats::rnorm(8, 0, 0.6)[dd$g]",
"  dd$pos <- exp(0.2 + 0.1 * dd$x + 0.3 * u + stats::rnorm(80, 0, 0.2))",
"  fit <- suppressWarnings(frm(bf(pos ~ x + (1 | g)), family = lognormal(),",
"                              data = dd))",
"  p <- lap_pair(fit)",
"  p$full$draws[3L, \"b_Intercept\"] <- 720",
"  p$lap$draws[3L, \"b_Intercept\"] <- 720",
"  a <- posterior_epred(p$lap, re_formula = NA)",
"  expect_true(any(is.infinite(a[3L, ])))",
"  expect_identical(a, posterior_epred(p$full, re_formula = NA))",
"})",
"")
writeLines(c(s[seq_len(i - 1L)], add, s[i:length(s)]), f)
cat("ok\n")
