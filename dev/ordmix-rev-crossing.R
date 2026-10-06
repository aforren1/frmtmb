# Reviewer of lane ordmix: brms 2.23.0's R-side category probabilities
# and posterior_predict draw at crossing thresholds (cs() offsets)
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
th <- matrix(c(0, 1, 0.5), 1)
p <- brms:::dcumulative(1:4, eta = 0, thres = th, disc = 1, link = "logit")
cat("dcumulative at thresholds 0, 1, 0.5:", format(p, digits = 4), "\n")
cp <- brms:::pordinal(1:4, eta = 0, disc = 1, thres = th,
                      family = "cumulative", link = "logit")
set.seed(1)
u <- runif(10000)
dr <- brms:::first_greater(matrix(cp, 10000, 4, byrow = TRUE), target = u)
cat("pordinal:", format(cp, digits = 4), "; posterior_predict's draw rule",
    "gives category shares", format(tabulate(dr, 4) / 1e4, digits = 3), "\n")
