# Lane ceplot punch 1: where a raw variable used only inside poly() or
# log(abs(z) + 1) lives in a fit, and what conditional_effects() does
# with it. Data seed 11.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(11)
d <- data.frame(x = rnorm(80), z = rnorm(80))
d$y <- rnorm(80, 1 + d$x - 0.5 * d$x^2 + log(abs(d$z) + 1))
fit <- frm(bf(y ~ poly(x, 2) + log(abs(z) + 1)), family = gaussian(), data = d)
cat("data_frame cols:", names(fit$frame$data_frame), "\n")
cat("fit$data cols:", names(fit$data), "\n")
r <- tryCatch(conditional_effects(fit), error = function(e) conditionMessage(e))
print(if (is.character(r)) r else names(r))
r <- tryCatch(conditional_effects(fit, "x"), error = function(e) conditionMessage(e))
print(if (is.character(r)) r else head(r$x[, c("x", "z", "estimate__")], 3))
r <- tryCatch(conditional_effects(fit, "z", data = d), error = function(e) conditionMessage(e))
print(if (is.character(r)) r else head(r$z[, c("x", "z", "estimate__")], 3))
