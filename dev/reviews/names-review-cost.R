# Review of lane names: cost of brms_variable_sort() on a fit with many
# latent states (ar(cov = FALSE), one zerr per row) and many r_ levels.
.libPaths(c("C:/Users/adf44/source/r/wt-names-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(9)
n <- 4000
d <- data.frame(x = rnorm(n), g = factor(rep(1:800, each = 5)),
                t = rep(1:5, 800))
d$y <- rpois(n, exp(0.3 + 0.2 * d$x + rnorm(800, 0, 0.4)[d$g]))
fit <- suppressWarnings(suppressMessages(
  frm(y ~ x + (1 + x | g) + ar(t, g), family = poisson(), data = d)))
nm <- brms_par_labels(fit, TRUE)
cat("names:", length(nm), "\n")
t1 <- system.time(for (i in 1:5) o <- brms_variable_sort(fit, nm))[[3]] / 5
cat(sprintf("brms_variable_sort: %.3f s per call on %d names\n", t1,
            length(nm)))
Rprof("dev/reviews/names-review-cost.prof", interval = 0.005)
for (i in 1:5) o <- brms_variable_sort(fit, nm)
Rprof(NULL)
print(head(summaryRprof("dev/reviews/names-review-cost.prof")$by.total, 12))
