# validate_prior(), get_prior(), brms prior frames, brmsfamily() and
# the family fields on the lane build. Seed 20260930.
# Output: dev/ordinal-log-smoke3.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
fam <- cumulative(threshold = "equidistant")
cat("fields:", fam$threshold, fam$link, fam$link_disc, "\n")
print(validate_prior(set_prior("normal(0, 1)", class = "delta"),
                     bf(y ~ x), family = fam, data = d))
print(get_prior(bf(y | thres(gr = g) ~ x), family = fam, data = d))
# brms's own prior frame, translated
bp <- brms::prior(normal(0, 1), class = delta) +
  brms::prior(normal(0, 3), class = Intercept)
fit <- frm(bf(y ~ x), family = fam, data = d, prior = bp)
print(fit$prior)
print(tryCatch(frm(bf(y | thres(gr = g) ~ x), family = fam, data = d,
                   prior = set_prior("normal(0, 1)", class = "delta",
                                     group = "c")),
               error = function(e) conditionMessage(e)))
fg <- frm(bf(y | thres(gr = g) ~ x), family = fam, data = d,
          prior = set_prior("normal(0.5, 0.01)", class = "delta",
                            group = "a"))
print(summary(fg)$spec_pars)
print(tryCatch(set_prior("normal(0, 1)", class = "delta", coef = "1") |>
                 (\(p) frm(y ~ x, family = fam, data = d, prior = p))(),
               error = function(e) conditionMessage(e)))
print(brmsfamily("acat", "probit", threshold = "sum_to_zero")$threshold)
print(tryCatch(cumulative(link_disc = "logit"), error = conditionMessage))
