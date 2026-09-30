# brms 2.23.0: what update.brmsformula() keeps when a complete bf()
# updates a stored one (the formula half of update.brmsfit()).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
old <- bf(count ~ 1 / (1 + exp(-a)) * exp(b * Trt), a ~ Age, b ~ Age,
          nl = TRUE)
print(update(old, bf(count ~ a + b, nl = TRUE), mode = "replace"))
print(update(old, bf(count ~ a + b, a ~ 1, nl = TRUE), mode = "replace"))
lin <- bf(y ~ x, sigma ~ z)
print(update(lin, bf(y ~ x2)))
print(update(lin, y ~ x2))
