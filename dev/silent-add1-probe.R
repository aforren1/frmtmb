# What add1() lists for a me() interaction on a fit that has the me().
.libPaths(c("C:/Users/adf44/source/r/wt-silent-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(11)
n <- 200
d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                k = sample(1:4, n, TRUE),
                g = factor(sample(letters[1:8], n, TRUE)))
d$xo <- d$x + stats::rnorm(n, 0, 0.3)
d$sdx <- 0.3
d$y <- 1 + 0.6 * d$k + 0.4 * d$z + sin(2 * d$x) +
  stats::rnorm(8)[d$g] + stats::rnorm(n, 0, 0.5)
fm <- frm(y ~ me(xo, sdx) + z, data = d)
print(add1(fm, scope = ~ . + me(xo, sdx):z, test = "Chisq"))
fi <- tryCatch(frm(y ~ me(xo, sdx) + z + me(xo, sdx):z, data = d),
               error = function(e) conditionMessage(e))
print(if (is.character(fi)) fi else logLik(fi))
print(logLik(fm))
