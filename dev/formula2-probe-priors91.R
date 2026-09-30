# priors:88-94's own case, with brms's answer beside it
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)),
                  x = rnorm(10), g = rep(1:2, 5))
family <- list(gaussian, Beta())
r <- tryCatch(default_prior(bf(y1 ~ x + (x|ID1|g) + ar()) + bf(y2 ~ 1),
                            dat, family = family),
              error = function(e) conditionMessage(e))
print(r)
pb <- brms::default_prior(brms::bf(y1 ~ x + (x|ID1|g) + ar()) +
                            brms::bf(y2 ~ 1), dat,
                          family = list(gaussian, brms::Beta()))
print(pb[, c("class", "coef", "group", "resp", "dpar")])
# the same model without ar(): the family list itself
dat2 <- dat; dat2$y2 <- dat$y2 / 4
r2 <- default_prior(bf(y1 ~ x + (x|ID1|g) + ar()) + bf(y2 ~ 1), dat2,
                    family = list(gaussian, Beta()))
print(r2)
