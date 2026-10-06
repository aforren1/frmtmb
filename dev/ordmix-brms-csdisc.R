# brms 2.23.0 on cs() in the formula of disc, and of mu2 of a mixture
# whose second component is cumulative. Output:
# dev/ordmix-log-brms-csdisc.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(y = sample(1:4, 60, TRUE), x = rnorm(60), z = rnorm(60))
r <- tryCatch(stancode(bf(y ~ x, disc ~ cs(z)), data = d, family = sratio()),
              error = function(e) e)
cat("disc ~ cs(z):", if (inherits(r, "error")) conditionMessage(r) else
  grep("cs", strsplit(as.character(r), "\n")[[1]], value = TRUE), sep = "\n")
r <- tryCatch(stancode(bf(y ~ x, mu1 ~ cs(z)), data = d,
                       family = mixture(sratio, cumulative)),
              error = function(e) e)
cat("mixture(sratio, cumulative), mu1 ~ cs(z):",
    if (inherits(r, "error")) conditionMessage(r) else
  grep("cs", strsplit(as.character(r), "\n")[[1]], value = TRUE), sep = "\n")
r <- tryCatch(stancode(bf(y ~ cs(z)), data = d,
                       family = mixture(sratio, cumulative)),
              error = function(e) e)
cat("mixture(sratio, cumulative), y ~ cs(z):",
    if (inherits(r, "error")) conditionMessage(r) else
  grep("cs", strsplit(as.character(r), "\n")[[1]], value = TRUE), sep = "\n")
