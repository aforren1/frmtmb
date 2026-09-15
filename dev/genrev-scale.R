LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...))
ulp <- function(a, b) {
  # distance in units in the last place, which is what "identical up to
  # rounding" means; %.6f renders 1.8e-11 as 0.000000
  max(abs(a - b) / pmax(.Machine$double.eps * abs(b), .Machine$double.xmin))
}
set.seed(2026)
n <- 400L
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))
f <- frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd)

mu <- predict(f, type = "link")
sg <- sigma(f)
ft <- fitted(f)
say("== the identity, at full precision ==\n")
say("max |fitted - exp(mu + s^2/2)| in ulps: %.3f\n", ulp(ft, exp(mu + sg^2 / 2)))
say("identical(fitted, exp(mu+s^2/2)):        %s\n",
    identical(ft, exp(mu + sg^2 / 2)))
r1 <- median(ft / exp(mu)); r2 <- exp(sg^2 / 2)
say("median(fitted/exp(mu)) = %.15f\n", r1)
say("exp(sigma^2/2)         = %.15f\n", r2)
say("difference in ulps     = %.3f\n", ulp(r1, r2))
say("max over rows of |fitted/exp(mu) - exp(s^2/2)| ulps = %.3f\n",
    ulp(ft / exp(mu), rep(r2, n)))

say("\n== the yardstick the lane replaced ==\n")
pl <- predict(f, type = "link", se.fit = TRUE)
se_pred <- pl$se.fit
# separation of the lognormal MEAN from its MEDIAN on the link scale is
# exactly sigma^2/2
sep_link <- sg^2 / 2
say("mean-minus-median on the link scale (s^2/2): %.6f\n", sep_link)
say("predict() own se.fit: min %.6f median %.6f max %.6f\n",
    min(se_pred), median(se_pred), max(se_pred))
say("separation / median se.fit = %.4f   <-- the lane's 0.9\n",
    sep_link / median(se_pred))
say("separation / min se.fit    = %.4f\n", sep_link / min(se_pred))
se_ls <- summary(f)$coefficients$sigma[1, 2]
se_ratio <- exp(sg^2 / 2) * sg^2 * se_ls
say("se of log(sigma) from summary(): %.6f\n", se_ls)
say("ratio = %.6f, its se = %.6f, t = %.4f\n", r2, se_ratio,
    (r2 - 1) / se_ratio)
say("the t statistic reduces to 1/(2*se_log_sigma) = %.4f, i.e. ~sqrt(n/2) = %.4f\n",
    1 / (2 * se_ls), sqrt(n / 2))

say("\n== does the stated lognormal identity survive a distributional sigma? ==\n")
f2 <- frm(bf(y ~ x + (1 | g)) + lf(sigma ~ x) + lognormal(), data = dd)
mu2 <- predict(f2, type = "link")
s2 <- sigma(f2)
say("sigma(f2) length = %d, value(s) = %s\n", length(s2),
    paste(signif(head(s2, 3), 6), collapse = ","))
sv <- predict(f2, type = "response", dpar = "sigma")
say("predict(dpar='sigma', type='response') length = %d range %.4f..%.4f\n",
    length(sv), min(sv), max(sv))
say("max |fitted(f2) - exp(mu + sigma(f2)^2/2)| relative = %.3e\n",
    max(abs(fitted(f2) - exp(mu2 + s2^2 / 2)) / fitted(f2)))
say("max |fitted(f2) - exp(mu + sv^2/2)| in ulps = %.3f\n",
    ulp(fitted(f2), exp(mu2 + sv^2 / 2)))

say("\n== does it survive truncation? ==\n")
dd2 <- dd[dd$y > 2000, ]
f3 <- tryCatch(frm(bf(y | trunc(lb = 2000) ~ x) + lognormal(), data = dd2),
               error = function(e) e)
if (inherits(f3, "error")) {
  say("trunc fit failed: %s\n", conditionMessage(f3))
} else {
  m3 <- predict(f3, type = "link"); s3 <- sigma(f3)
  say("max |fitted - exp(mu+s^2/2)| relative = %.4e  (0 would mean the\n",
      max(abs(fitted(f3) - exp(m3 + s3^2 / 2)) / fitted(f3)))
  say("  page's lognormal formula still holds under truncation)\n")
}

say("\n== what the contract page does NOT state ==\n")
say("conditional_effects() estimate scale, row 1:\n")
ce <- tryCatch(conditional_effects(f, effects = "x"), error = function(e) e)
if (!inherits(ce, "error")) {
  d1 <- ce[[1]]
  say("  columns: %s\n", paste(names(d1), collapse = ","))
  say("  estimate range %.4f .. %.4f  (link ~8, response ~1e3)\n",
      min(d1$estimate__), max(d1$estimate__))
} else say("  error: %s\n", conditionMessage(ce))
ps <- posterior_summary(f)
say("posterior_summary(fit) Estimate[1:2] = %s\n",
    paste(signif(head(ps[, 1], 2), 8), collapse = ","))
say("bayes_R2 / hypothesis / pp_check scale: not stated on the page\n")
cat("GENREVDONE\n")
