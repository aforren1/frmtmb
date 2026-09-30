args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] == "before") character() else
  "C:/Users/adf44/source/r/wt-formula2-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("==", label, "\n")
  if (is.character(r)) cat(r, "\n") else print(r)
}
dat <- data.frame(y = 1:10, x1 = 1:10, x2 = (1:10)^2)
tr("dot frame", colnames(frm(y ~ ., data = dat, dry_run = "frame")$linpreds[[1]]$X))
tr("dot fit", fixef(frm(y ~ ., data = dat)))
tr("dot sigma", frm(bf(y ~ x1, sigma ~ .), data = dat, dry_run = "frame")$linpreds[[2]]$X[1:2, ])
# family list
d2 <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), x = rnorm(10),
                 g = rep(1:2, 5))
tr("family list", get_prior(bf(y1 ~ x) + bf(y2 ~ 1), data = d2,
                            family = list(gaussian, poisson())))
# cmc
df <- data.frame(y = 1:10, g = rep(c("a", "b"), 5))
tr("lf cmc", lf(disc ~ 0 + g, cmc = FALSE))
tr("bf cmc", bf(y ~ 0 + g, cmc = FALSE))
# equate
tr("equate", bf(y ~ x, sigma1 = "sigma2"))
