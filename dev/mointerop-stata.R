## What marginaleffects 1.0.0 gives lm() for vcov = "stata", "robust", "HC1"
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(marginaleffects))
set.seed(31); n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:20, n, TRUE)))
d$y <- 1 + 0.5 * d$x + rnorm(n) * exp(0.3 * d$x)
m <- lm(y ~ x, data = d)
for (v in c("HC0", "HC1", "stata", "robust", "HC3")) {
  cat(v, sqrt(diag(get_vcov(m, vcov = v))), "\n")
}
cat("vcovHC HC1", sqrt(diag(sandwich::vcovHC(m, type = "HC1"))), "\n")
cat("vcovHC HC3", sqrt(diag(sandwich::vcovHC(m, type = "HC3"))), "\n")
gm <- glm(y ~ x, data = d)
for (v in c("HC0", "HC1", "stata")) cat("glm", v, sqrt(diag(get_vcov(gm, vcov = v))), "\n")
