# Lane optrap, punch round 1: the unit test-se-check.R's m3 data gets
# under the residual-spread rule, as drawn and halved.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(101)
d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
for (k in c(1, 2)) {
  d <- d3
  d$y <- d$y / k
  fr <- suppressWarnings(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
                             data = d, prior = set_prior("", class = "b",
                                                         ub = 0.1)))
  cat("y /", k, "units:", deparse(frmtmb:::predictor_units(fr$frame)),
      "sd(resid lm)", sd(residuals(lm(y ~ x, d))), "\n")
}
