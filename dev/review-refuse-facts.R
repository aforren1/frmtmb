# Small facts the review quotes. Run twice: REV_LIB unset (lane) and with
# REV_LIB=C:/Users/adf44/source/r/rellib-r16 (before). Seed 10800.
source("dev/review-refuse-lib.R")
pr("build: %s", dirname(find.package("frmtmb")))
set.seed(10800); n <- 120; x <- rnorm(n); g <- factor(sample(1:12, n, TRUE))
d <- data.frame(x = x, g = g, y = rnorm(n, 0.5 * x))
a <- tr(frm(bf(y ~ x + (1 | gr(g, id = "p")), sigma ~ (1 | gr(g, id = "p", cor = FALSE))), data = d))
pr("id mixed flags: %s", msg(a, 160))
a <- tr(frm(bf(y ~ x + (1 | gr(g, id = "p", cor = FALSE))), data = d))
pr("id cor = FALSE alone: %s", msg(a, 160))
# cooks.distance at the missing rows
mrows <- seq(4, n, 12)
d$y[mrows] <- NA
f <- frm(bf(y | mi() ~ x), data = d)
ck <- cooks.distance(f)
pr("gaussian mi(): cooks.distance at missing rows %s; at others, max %.4f", paste(signif(ck[mrows][1:4], 3), collapse = " "), max(ck[-mrows], na.rm = TRUE))
cc <- frm(bf(y ~ x), data = d[-mrows, ])
pr("cooks of the complete-case fit, max %.4f", max(cooks.distance(cc), na.rm = TRUE))
if (dirname(find.package("frmtmb")) != "C:/Users/adf44/source/r/rellib-r16") {
  d2 <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))); d2$y[mrows] <- NA
  f2 <- frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d2)
  ck2 <- cooks.distance(f2)
  pr("Gamma mi(): cooks.distance at missing rows %s", paste(signif(ck2[mrows][1:4], 3), collapse = " "))
}
# the s(I(x^2)) refusal
d3 <- data.frame(x = runif(n), y = rnorm(n))
a <- tr(frm(bf(y ~ s(I(x^2), k = 5)), data = d3))
pr("s(I(x^2)): %s", msg(a, 120))
