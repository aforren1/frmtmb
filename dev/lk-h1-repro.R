suppressMessages(library(frmtmb))
options(digits = 8)
set.seed(11); n <- 400; x <- rnorm(n)
d <- data.frame(x = x)
mig <- 1 / sqrt(pmax(0.5 + 0.2 * x, 0.05))
d$yig <- rgamma(n, shape = 20, rate = 20 / mig)
f <- frm(yig ~ x, data = d, family = inverse.gaussian(link = "1/mu^2"))
z <- conditional_effects(f)[[1]]
cat("inverse.gaussian 1/mu^2\n")
cat("  inverted rows:", sum(z$lower__ > z$upper__, na.rm = TRUE), "of", nrow(z), "\n")
cat("  non-finite bound rows:",
    sum(!is.finite(z$lower__) | !is.finite(z$upper__)), "\n")
print(utils::head(data.frame(est = z$estimate__, lo = z$lower__,
                             up = z$upper__), 3))

cat("\nGamma(inverse), pre-existing on unchanged code\n")
d$yg <- rgamma(n, shape = 20, rate = 20 / (1 / pmax(0.5 + 0.2 * x, 0.05)))
fg <- frm(yg ~ x, data = d, family = Gamma(link = "inverse"))
zg <- conditional_effects(fg)[[1]]
cat("  inverted rows:", sum(zg$lower__ > zg$upper__, na.rm = TRUE), "of",
    nrow(zg), "\n")
cat("  non-finite bound rows:",
    sum(!is.finite(zg$lower__) | !is.finite(zg$upper__)), "\n")

cat("\ncontrols\n")
d$yl <- rpois(n, exp(0.5 + 0.3 * x))
fl <- frm(yl ~ x, data = d, family = poisson())
zl <- conditional_effects(fl)[[1]]
cat("  poisson(log) inverted:", sum(zl$lower__ > zl$upper__, na.rm = TRUE), "\n")

cat("\npredict(se.fit = TRUE) response-scale interval, 1/mu^2\n")
nd <- data.frame(x = c(-2, -1, 0, 1, 2))
pr <- stats::predict(f, newdata = nd, type = "response", se.fit = TRUE)
pl <- stats::predict(f, newdata = nd, type = "link", se.fit = TRUE)
print(data.frame(eta = pl$fit, se_eta = pl$se.fit,
                 mu = pr$fit, se_mu = pr$se.fit))
