# Reviewer: printing of equations and the Links line with two sigma links
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
print(bf(y ~ x, sigma1 = "sigma2"))
print(lf(sigma1 = "sigma2", sigma ~ z))
print(bf(y ~ 0 + g, cmc = FALSE))
set.seed(11)
d <- data.frame(x = rnorm(300))
k <- rbinom(300, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(300, 3 + 0.5 * d$x, 1), rnorm(300, -1, 1))
f3 <- suppressMessages(frm(bf(y ~ x, sigma1 = "sigma2"), data = d,
  family = mixture(brmsfamily("gaussian", link_sigma = "softplus"),
                   gaussian())))
s <- capture.output(print(f3))
cat(s[1:4], sep = "\n")
print(frm(bf(y ~ x, sigma1 = "sigma2"), data = d,
          family = mixture(gaussian(), gaussian()), dry_run = "spec"))
