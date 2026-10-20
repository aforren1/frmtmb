# brms 2.23.0's update() refusal of a variable its stored data lack,
# word for word, from an empty fit with a stub stanfit (no compile; the
# stub is the reviewer's construction in dev/drop1raw-rev-brms2.R).
# Log: dev/drop1raw-out/brms-newvar.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages({library(brms); library(rstan)})
cat("brms", format(packageVersion("brms")), "\n")
set.seed(5)
d <- data.frame(x = rnorm(40), x2 = rexp(40) + 0.1, w = rnorm(40),
                v = rnorm(40))
d$z <- 1 + d$x + rnorm(40)
b <- brm(z ~ x + log(x2), data = d, empty = TRUE)
b$fit <- methods::new("stanfit", sim = list(
  iter = 2, warmup = 1, chains = 1, thin = 1,
  samples = list(structure(list(), args = list(control = list())))))
for (f in list(~ . + w, ~ . + w + v)) {
  e <- tryCatch(update(b, f, recompile = TRUE, empty = TRUE),
                error = function(e) e)
  cat(deparse(f), "->\n")
  print(conditionMessage(e))
}
