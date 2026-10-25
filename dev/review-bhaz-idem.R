# Reviewer: cox_bhaz_option() applied twice, and to a cox() that set df.
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.survival))
op <- frmtmb.survival:::cox_bhaz_option
set.seed(2); n <- 300; x <- rnorm(n)
tt <- stats::rexp(n, exp(-0.5 + 0.7 * x)); ct <- stats::rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
for (cl in list(quote(bhaz(df = 6)), quote(bhaz(knots = c(0.5, 1.1))),
                quote(bhaz(degree = 2, intercept = FALSE)))) {
  f1 <- op(cox(), cl, globalenv())
  f2 <- op(f1, cl, globalenv())
  fit <- function(fam) frm(bf(time | cens(cens) ~ x), family = fam, data = dd)
  a <- fit(f1); b <- fit(f2)
  cat(deparse1(cl), ": cox_args identical", identical(f1$cox_args, f2$cox_args),
      " fits identical", identical(logLik(a), logLik(b)) &&
        identical(a$estimates, b$estimates), "\n")
}
