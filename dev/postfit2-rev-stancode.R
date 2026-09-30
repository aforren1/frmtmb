.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
set.seed(22)
db <- data.frame(x = rnorm(160), g = factor(rep(1:8, 20)), h = factor(rep(1:10, each = 16)), y = rnorm(160))
sc <- strsplit(stancode(y ~ x + (1 | g) + (1 | h), data = db), "\n")[[1]]
cat(paste(seq_along(sc), sc)[25:45], sep = "\n")
