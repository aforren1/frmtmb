# Reviewer probe 6: brms 2.23.0 on b/Intercept spellings with
# center = FALSE and a dpar formula. Seed 1.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(x = rnorm(40)); d$y <- rnorm(40)
msg <- function(expr) tryCatch({ force(expr); "ok" },
                               error = function(e) conditionMessage(e))
f <- bf(y ~ x, sigma ~ x, center = FALSE)
cat("brms validate b/Intercept dpar sigma, center = FALSE:",
    msg(validate_prior(set_prior("normal(0, 1)", class = "b",
                                 coef = "Intercept", dpar = "sigma"),
                       f, data = d)), "\n")
cat("brms validate Intercept dpar sigma, center = FALSE:",
    msg(validate_prior(set_prior("normal(0, 1)", class = "Intercept",
                                 dpar = "sigma"), f, data = d)), "\n")
tab <- default_prior(f, data = d)
print(tab[, c("prior", "class", "coef", "dpar")])
sc <- stancode(f, data = d)
cat(grep("Intercept", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
cat("\nbrms: lb on class sd:", msg(validate_prior(
  set_prior("", class = "sd", lb = 0.1), bf(y ~ x + (1 | g)),
  data = transform(d, g = rep(1:8, 5)))), "\n")
