# Reviewer of lane fixes, claim 4b: what brms 2.23.0's Intercept is
# under threshold = "sum_to_zero" (stancode only, no compile).
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(1)
n <- 100
d <- data.frame(x = rnorm(n))
d$y <- 1L + (d$x + rlogis(n) > -1) + (d$x + rlogis(n) > 0) +
  (d$x + rlogis(n) > 1)
for (fam in c("sratio", "cumulative")) {
  sc <- brms::stancode(y ~ x, data = d,
                       family = get(fam, asNamespace("brms"))(
                         threshold = "sum_to_zero"),
                       prior = brms::prior(normal(0, 1), class = Intercept,
                                           coef = "2"))
  cat("==", fam, "\n")
  l <- strsplit(sc, "\n")[[1]]
  print(grep("Intercept|stz|sum_to_zero", l, value = TRUE))
}
