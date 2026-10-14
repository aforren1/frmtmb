# Reviewer of lane illse: parameter names and maps (exploration only).
.libPaths(c("C:/Users/adf44/source/r/wt-illse-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(1)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, length.out = 200)))
d$y <- sin(3 * d$x) + rnorm(20, 0, .5)[d$g] + rnorm(200, 0, .3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
f <- frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d)
print(f$opt$par)
print(names(f$obj$env$par))
print(f$opt$convergence)
print(str(ns$sdr_of(f)$cov_fixed_root))
print(names(f$estimates))
# cumulative
d$yo <- factor(cut(2 * sin(3 * d$x) + rlogis(200), c(-Inf, -1, 0, 1, Inf),
                   labels = FALSE), ordered = TRUE)
fc <- frm(yo ~ x1 + x2 + x3 + x4 + x5 + x6 + x7, data = d,
          family = cumulative())
print(fc$opt$par)
post <- fc$spec$responses[[1]]$family$post
print(post$ord_thresholds)
print(post$ord_thresholds_raw)
print(fixef(fc))
print(args(ns$predict.frmtmb_fit))
print(args(ns$frm_multiple))
print(args(ns$confint.frmtmb_fit))
