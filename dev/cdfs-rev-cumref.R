# Reviewer: which refusal the fuzz probe cens_no_cdf / trunc_no_cdf now
# meets with cumulative() (it needs "CDF" in the message)
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(60))
d$y <- factor(sample(1:4, 60, TRUE), ordered = TRUE)
d$cc <- sample(c(-1, 0, 1), 60, TRUE)
for (f in list(y | cens(cc) ~ x, y | trunc(lb = 2) ~ x)) {
  r <- tryCatch({frm(f, family = cumulative(), data = d); "FITTED"},
                error = function(e) conditionMessage(e))
  cat(deparse(f), "->", r, "\n\n")
}
