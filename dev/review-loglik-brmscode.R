# Review of lane loglik: brms 2.23.0's lprior lines for a half prior on
# sd and a lower-bounded prior on b (stancode only, no compile).
.libPaths(c("C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("lib: brms", format(packageVersion("brms")), dirname(find.package("brms")),
    "\n")
set.seed(304)
d <- data.frame(x = rnorm(60), g = factor(rep(1:12, each = 5)))
d$y <- 1 + 0.5 * d$x + rnorm(60)
for (p in list(prior(normal(0, 1), class = sd),
               prior(normal(0, 1), class = b, lb = 0))) {
  s <- strsplit(stancode(y ~ x + (1 | g), data = d, prior = p), "\n")[[1]]
  i <- grep("lprior +=", s, fixed = TRUE)
  cat(s[sort(unique(c(i, i + 1)))], sep = "\n")
  cat("--\n")
}
