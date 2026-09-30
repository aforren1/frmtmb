# Lane wt-defects, punch round 1, B2: what brms 2.23.0's
# fitted(scale = "linear") returns on a multivariate fit with an ordinal
# cs() response, on the reviewer's construction (dev/defects-rev-probe5.R,
# seed 5, n = 150), and on the ordinal response alone. A short chain:
# only the shapes and names are read. Log dev/defects-log/brms-mvcs.txt.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
set.seed(5)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$yo <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
d$yg <- d$z + rnorm(n)
f <- brm(bf(yg ~ x) + bf(yo ~ x + cs(z), family = sratio()) +
           set_rescor(FALSE), data = d, chains = 1, iter = 300,
         refresh = 0, seed = 1)
a <- fitted(f, scale = "linear")
cat("multivariate dim:", dim(a), "\n")
print(dimnames(a))
u <- brm(yo ~ x + cs(z), family = sratio(), data = d, chains = 1,
         iter = 300, refresh = 0, seed = 1)
b <- fitted(u, scale = "linear")
cat("univariate dim:", dim(b), "\n")
print(dimnames(b))
