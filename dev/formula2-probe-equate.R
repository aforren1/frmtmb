.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.4)
d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
f <- bf(y ~ x, sigma1 = "sigma2")
str(unclass(f))
fam <- mixture(gaussian, gaussian)
cat(stancode(f, data = d, family = fam))
print(default_prior(f, data = d, family = fam))
# other spellings
f2 <- bf(y ~ x) + lf(sigma1 = "sigma2")
print(try(str(unclass(f2))))
f3 <- try(bf(y ~ x, flist = list(sigma1 = "sigma2")))
print(f3)
# theta equated? theta class
print(try(bf(y ~ x, theta1 = "theta2")), silent = TRUE)
# theta equating crashes brms in bf()
# mu class refused; with itself refused
print(try(bf(y ~ x, mu1 = "mu2")))
# fixed plus equated: sigma2 = 1 and sigma1 = "sigma2"
print(try(bf(y ~ x, sigma2 = 1, sigma1 = "sigma2")))
# bf updated: bf(bf(y ~ x, sigma1 = "sigma2"), sigma2 ~ x)
print(try(bf(bf(y ~ x, sigma1 = "sigma2"), sigma2 ~ x)))
# non-mixture: student nu? Equate among different classes refused
print(try(bf(y ~ x, sigma = "nu")))
# brmsterms on something where rhs not a valid dpar
print(try(brmsterms(bf(y ~ x, sigma1 = "sigma3"), family = fam)))
print(try(stancode(bf(y ~ x, sigma1 = "sigma3"), data = d, family = fam)))
# vector-valued character
print(try(bf(y ~ x, sigma1 = c("sigma2", "sigma3"))))
# a character for a non-mixture family
print(try(stancode(bf(y ~ x, sigma = "sigma"), data = d)))
print(try(bf(y ~ x, shape1 ~ x, shape2 = "shape1")))
print(try(bf(y ~ x, shape1 = "shape3", shape2 = "shape1")))
# three components: sigma1 = sigma2 = sigma3 via two equations
f4 <- bf(y ~ x, sigma1 = "sigma3", sigma2 = "sigma3")
cat(stancode(f4, data = d, family = mixture(gaussian, gaussian, gaussian)))
