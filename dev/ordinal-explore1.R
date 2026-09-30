# Exploration on the BASE build (rellib-r4): how hurdle_cumulative's disc
# and the ordinal thresholds surface in the post-fit methods today.
.libPaths(c("C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
set.seed(18)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
yc <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
d$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * d$x), 0L, yc)
d$yo <- yc
fit <- frm(bf(y ~ x, disc ~ 0 + z), family = hurdle_cumulative(), data = d)
print(fit)
cat("-- variables --\n"); print(variables(fit))
cat("-- fixef --\n"); print(fixef(fit))
cat("-- summary --\n"); print(summary(fit))
cat("-- default_prior --\n")
print(default_prior(bf(y ~ x, disc ~ z), family = hurdle_cumulative(),
                    data = d))
fit0 <- frm(bf(yo ~ x), family = cumulative(), data = d)
cat("-- cumulative variables --\n"); print(variables(fit0))
print(summary(fit0))
cat("-- fitted linear --\n")
print(head(fitted(fit0, scale = "linear")))
print(dim(fitted(fit0)))
print(tryCatch(frm(bf(yo ~ x, disc ~ z), family = cumulative(), data = d),
               error = function(e) conditionMessage(e)))
print(tryCatch(cumulative(threshold = "equidistant"),
               error = function(e) conditionMessage(e)))
print(args(cumulative)); print(args(sratio))
