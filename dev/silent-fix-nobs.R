# Fix round, MINOR 4: which count each reader of n gives on a
# multivariate subset() fit (rows 1..72 of 80), on the lane build.
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5)
d6 <- data.frame(x = rnorm(80))
d6$y1 <- d6$x + rnorm(80)
d6$y2 <- d6$x + rnorm(80)
d6$s1 <- seq_len(80) <= 72
f <- frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s1) ~ x) +
           set_rescor(FALSE), data = d6)
ll <- logLik(f)
cat("nobs()", nobs(f), " insight::n_obs()",
    tryCatch(insight::n_obs(f), error = function(e) conditionMessage(e)),
    " logLik nobs", attr(ll, "nobs"),
    " BIC - AIC(k = log(72))",
    format(BIC(f) - AIC(f, k = log(72)), digits = 3), "\n")
