# conditional_effects() on gp(x) and gp(x, by = f): which effects.
lib <- Sys.getenv("GPBY_LIB", "C:/Users/adf44/source/r/wt-gpby-lib")
.libPaths(c(if (nzchar(lib) && lib != "base") lib,
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 90
dd <- data.frame(x = runif(n, 0, 5),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)))
dd$y <- ifelse(dd$f == "a", sin(dd$x), cos(dd$x)) + rnorm(n, 0, 0.3)
f0 <- frm(bf(y ~ gp(x)), data = dd)
ce0 <- conditional_effects(f0)
print(names(ce0))
str(lapply(ce0, head, 2))
f1 <- try(frm(bf(y ~ gp(x, by = f)), data = dd))
if (!inherits(f1, "try-error")) {
  ce1 <- conditional_effects(f1)
  print(names(ce1))
  print(head(ce1[["x:f"]]))
}
