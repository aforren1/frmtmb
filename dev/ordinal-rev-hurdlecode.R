# Reviewer: brms's generated hurdle_cumulative code with a modeled disc.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
set.seed(1); n <- 50
d <- data.frame(x = rnorm(n), z = rnorm(n), y = sample(0:5, n, TRUE))
sc <- brms::make_stancode(brms::bf(y ~ x, disc ~ 0 + z), data = d,
  family = brms::hurdle_cumulative(threshold = "equidistant"))
L <- strsplit(sc, "\n")[[1]]
cat(paste(seq_along(L), L)[c(10:40, 100:120)], sep = "\n")
