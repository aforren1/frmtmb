# Lane sampfix nits round: a test that reaches ce_at() with new levels.
f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/extensions/frmtmb.sample/tests/testthat/test-laplace-draws.R"
s <- gsub("\r", "", readLines(f))
i <- which(s == "test_that(\"mi() values are integrated too, and only their reader refuses\", {")
stopifnot(length(i) == 1L)
add <- c(
"test_that(\"a curve at re_formula = NULL that reads a smooth is refused\", {",
"  # re_formula = NULL draws a new group level per draw, and the smooth's",
"  # coefficients are still read: the refusal must see that, not an",
"  # error of its own that it takes for one unrelated to the fill (R CMD",
"  # check found such a probe calling a helper that no longer existed,",
"  # and the curves then came back NA without a word)",
"  dd <- lap_data()",
"  fit <- suppressWarnings(frm(bf(y ~ s(x, k = 5) + (1 | g)),",
"                              family = gaussian(), data = dd))",
"  p <- lap_pair(fit)",
"  expect_error(conditional_effects(p$lap, effects = \"x\", resolution = 5,",
"                                   re_formula = NULL, seed = 2),",
"               \"conditional_effects() needs the random effects\",",
"               fixed = TRUE)",
"})",
"")
writeLines(c(s[seq_len(i - 1L)], add, s[i:length(s)]), f)
cat("ok\n")
