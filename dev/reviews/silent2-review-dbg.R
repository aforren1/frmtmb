# Review of lane silent2: diagnose factor(cyl) contrasts and I(y) anova.
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib", "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
set.seed(7)
dc <- data.frame(f = factor(sample(c("a", "b", "c"), 90, TRUE)),
                 g = factor(rep(1:10, each = 9)), x = rnorm(90),
                 cyl = sample(c(4, 6, 8), 90, TRUE))
dc$y <- rnorm(10)[dc$g] + rnorm(90)
fc <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                          q(y ~ factor(cyl) + (1 + factor(cyl) | g), data = dc))
for (bk in fc$frame$re_blocks) for (cp in bk$components) {
  cat("component contrasts:\n"); str(cp$contrasts)
}
cat("cnms:", fc$frame$re_blocks[[1]]$cnms, "\n")
print(tryCatch(fitted(fc, newdata = dc[1:3, ]), error = function(e) conditionMessage(e)))
# what re_term_contrasts sees
print(frmtmb:::re_term_contrasts(quote(1 + factor(cyl)), dc, globalenv()))
# anova I(y)
d <- data.frame(x = rnorm(50)); d$y <- d$x + rnorm(50)
a <- q(y ~ 1, data = d); b <- q(I(y) ~ x, data = d)
str(frmtmb:::loglik_rows_key(a)); str(frmtmb:::loglik_rows_key(b))
b2 <- q(y ~ x, data = d)
cat("y~1 vs y~x:", identical(frmtmb:::loglik_rows_key(a), frmtmb:::loglik_rows_key(b2)), "\n")
# two univariate fits of different responses with the same values
d$y2 <- d$y
c2 <- q(y2 ~ x, data = d)
print(tryCatch(anova(a, c2), error = function(e) conditionMessage(e)))
