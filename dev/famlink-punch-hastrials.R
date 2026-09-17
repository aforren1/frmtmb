# Which brms families require trials(), by brms's own predicate.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
print(brms:::has_trials)
fs <- sub("^[.]family_", "", ls(asNamespace("brms"), all.names = TRUE,
                                pattern = "^[.]family_"))
fs <- setdiff(fs, c("custom", "info"))
cat("has_trials:", fs[vapply(fs, function(f) isTRUE(brms:::has_trials(f)), NA)], "\n")
d <- data.frame(y = rbinom(20, 1, 0.5), x = rnorm(20))
r <- tryCatch(brms::make_standata(y ~ x, d, family = brms::mixture(binomial, binomial)),
              error = conditionMessage)
cat("mixture(binomial, binomial) without trials:", if (is.character(r)) r else "OK", "\n")
