# Reviewer: how brms 2.23.0 groups NEW-data gp() rows. The lane ran
# standata() on fitting data only. Here an empty brmsfit (no compile)
# and standata(fit, newdata = ) on unseen rows one ulp apart.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(17)
xg <- round(runif(60, 0, 10), 1)
d <- data.frame(y = sin(xg) + rnorm(60, 0, 0.3), x = xg)
fit <- brm(y ~ gp(x), data = d, empty = TRUE)
u <- 10 / 3
nd <- data.frame(x = c(u, u * (1 + 2^-52), u, 1 / 3, 1 / 3 * (1 + 2^-52)),
                 y = 0)
cat("distinct doubles in newdata x:", length(unique(nd$x)), "\n")
sd <- standata(fit, newdata = nd, internal = TRUE)
cat("Jgp_1 (new rows -> new positions):", sd$Jgp_1, "\n")
cat("Nsubgp_1:", sd$Nsubgp_1, " rows of Xgp_1:", NROW(sd$Xgp_1), "\n")
cat("rows of Xgp_1_old (fitted positions):", NROW(sd$Xgp_1_old), "\n")
