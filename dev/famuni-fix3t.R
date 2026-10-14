p <- "tests/testthat/test-families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
old <- '    allowed = c("did not report convergence", "maximum absolute gradient",
                "Standard errors"))'
new <- '    allowed = c("did not report convergence", "maximum absolute gradient",
                "Standard errors", "the edge at -1"),
    require = "the edge at -1")'
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
x <- paste0(x, '

test_that("the xi edge warning does not fire inside the range", {
  # the case where the guard must stay quiet: xi = -0.45 and 0.8 at
  # n = 150 end at -0.478 and 0.889 (the review\'s section A)
  rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
  for (cfg in list(c(-0.45, 401), c(0.8, 406))) {
    set.seed(cfg[2])
    d <- data.frame(y = rg(150, 2, 1.5, cfg[1]))
    expect_no_warning(fit <- frm(y ~ 1, family = gen_extreme_value(),
                                 data = d))
    xi <- frm_linpred(fit, dpar = "xi", type = "response")[1]
    expect_gt(xi, -0.9)
  }
})
')
writeLines(strsplit(x, "\n")[[1]], p)
