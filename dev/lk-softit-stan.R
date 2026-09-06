set.seed(1)
d <- data.frame(y = rbinom(30, 1, 0.5), x = rnorm(30))
sc <- as.character(brms::stancode(brms::bf(y ~ x), data = d,
                                  family = brms::bernoulli(link = "softit")))
ln <- strsplit(sc, "\n")[[1]]
i <- grep("softit", ln)
cat(ln[max(1, min(i) - 3):(max(i) + 1)], sep = "\n")
cat("\n--- the offending line, exactly:\n")
bad <- grep("expm1", ln, value = TRUE)
print(bad)
cat("\n--- does the vector form use / where Stan needs ./ ?\n")
print(any(grepl("expm1(-p / (p - 1))", ln, fixed = TRUE)))
