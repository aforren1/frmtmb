# Reviewer: (cs(1) | g) with and without brms attached, rellib-r6 and r5.
for (lib in c("C:/Users/adf44/source/r/rellib-r6", "C:/Users/adf44/source/r/rellib-r5")) {
  cat(lib, "\n")
  cmd <- sprintf('.libPaths(c("%s", "C:/Users/adf44/AppData/Local/R/win-library/4.6")); %s suppressPackageStartupMessages(library(frmtmb)); set.seed(1); d <- data.frame(x = rnorm(200), g = factor(sample(letters[1:6], 200, TRUE))); d$y <- sample(1:4, 200, TRUE); cat(tryCatch({frm(y ~ x + (cs(1) | g), family = sratio(), data = d); "FIT"}, error = function(e) conditionMessage(e)), "\n")', lib, "%s")
  for (pre in c("", "suppressPackageStartupMessages(library(brms));")) {
    f <- tempfile(fileext = ".R"); writeLines(sprintf(cmd, pre), f)
    cat("  brms attached:", nzchar(pre), ":", system2(file.path(R.home("bin"), "Rscript.exe"), f, stdout = TRUE, stderr = TRUE), "\n")
  }
}
