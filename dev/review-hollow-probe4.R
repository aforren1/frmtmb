# Reviewer probe 4: make_conditions() on draws on the BASE build
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(2)
dd <- data.frame(x = rnorm(40))
dd$y <- 1 + 0.5 * dd$x + rnorm(40)
fp <- suppressWarnings(frm(y ~ poly(x, 2), data = dd))
fl <- suppressWarnings(frm(y ~ x, data = dd))
for (f in list(fp, fl)) {
  dp <- structure(list(stanfit = NULL, draws = matrix(0, 1, 1), fit = f),
                  class = "frmtmb_draws")
  r <- tryCatch(make_conditions(dp, "x"), error = function(e) e)
  cat("base make_conditions(draws, 'x'):",
      if (inherits(r, "error")) conditionMessage(r) else
        paste(format(unlist(r$x), digits = 4), collapse = " "), "\n")
}
