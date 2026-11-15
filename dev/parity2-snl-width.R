# Lane parity2, fix round F1/F2: the fit's band widths at a new level on
# the review's strong-shrinkage data (30 groups of 2 rows, sd_g 1,
# sigma 2, data seed 100), x = -1, 0, 1, re_formula = NULL.
#   Rscript dev/parity2-snl-width.R > dev/parity2-out/snl-width.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(100)
d <- data.frame(x = rnorm(60), g = factor(rep(1:30, each = 2)))
d$y <- 1 + 0.5 * d$x + rnorm(30)[d$g] + rnorm(60, sd = 2)
fb <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
ic <- list(x = c(-1, 0, 1))
wid <- function(ce) ce$x$upper__ - ce$x$lower__
for (B in c(200, 400)) {
  cb <- suppressMessages(conditional_effects(fb, "x", re_formula = NULL,
                                             band = "boot", boot = B,
                                             seed = 2, int_conditions = ic))
  cat("default boot", B, "seed 2 width:", format(wid(cb), digits = 4), "\n")
}
cw <- conditional_effects(fb, "x", re_formula = NULL, int_conditions = ic)
cat("wald width:", format(wid(cw), digits = 4), "\n")
for (s in c("uncertainty", "old_levels")) {
  for (b in c("boot", "wald")) {
    e <- tryCatch(conditional_effects(fb, "x", re_formula = NULL, band = b,
                                      boot = 4, int_conditions = ic,
                                      sample_new_levels = s),
                  error = function(e) conditionMessage(e))
    cat(s, b, ":", substr(if (is.character(e)) e else "answered", 1, 120),
        "\n")
  }
}
