# brms 2.23.0: c() of a set_prior() and a default_prior() table, in both
# orders; what the result is and what brm() reads from it (stancode()).
# Seed 11, the release review's y ~ x + (1 + x | g).
#   Rscript dev/rel077-fix2-brms-c.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(11)
d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] + rnorm(60)
f <- y ~ x + (1 + x | g)
bp <- suppressMessages(default_prior(f, data = d))
sp <- set_prior("normal(0, 1)", class = "b")
for (nm in c("c(sp, bp)", "c(bp, sp)")) {
  r <- tryCatch(eval(parse(text = nm)), error = function(e) e)
  if (inherits(r, "error")) { cat(nm, "ERROR:", conditionMessage(r), "\n"); next }
  cat(nm, "class", class(r)[1], "rows", nrow(r), "\n")
  print(as.data.frame(r)[, c("prior", "class", "coef", "group", "source")])
  sc <- tryCatch(stancode(f, data = d, prior = r), error = function(e) e)
  if (inherits(sc, "error")) cat("  stancode ERROR:", conditionMessage(sc), "\n") else
    cat("  stancode lprior lines:", grep("lprior", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n    ")
}
