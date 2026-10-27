## brms 2.23.0 on the fix-round-1 cases of dev/mocode-fix1.R (same data):
## conditional_effects() reference and moderator values, emmeans, the
## posterior sd of fitted values (the simplex's share of a prediction's
## uncertainty), and the empty-category case. 2 chains of 2000, seed 3.
##   Rscript dev/mocode-brms-fix1.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + 0.5 * d$z + rnorm(n)
fz <- q(brm(y ~ mo(k) + z, data = d, chains = 2, iter = 2000, seed = 3,
            refresh = 0))
ce <- q(conditional_effects(fz, effects = "z"))[[1]]
cat("effects = z: k held at", unique(ce$k), "; first z", ce$z[1],
    "estimate", signif(ce$estimate__[1], 7), "\n")
ce2 <- q(conditional_effects(fz, effects = "z:k"))[[1]]
cat("effects = z:k: k values", unique(ce2$k), "; effect2__ labels",
    unique(as.character(ce2$effect2__)), "\n")
ce3 <- q(conditional_effects(fz, effects = "k"))[[1]]
cat("effects = k: grid", ce3$k, "; z held at", unique(ce3$z), "\n")
r <- tryCatch(q(emmeans::emmeans(fz, ~ z)), error = function(e) e)
cat("emmeans(~ z):", if (inherits(r, "error")) paste("ERROR:",
    conditionMessage(r)) else format(summary(r)$emmean), "\n")
r <- q(emmeans::emmeans(fz, ~ z, at = list(k = 3)))
cat("emmeans(~ z, at = list(k = 3)):", format(summary(r)$emmean), "\n")
f <- fitted(fz, newdata = data.frame(k = 1:4, z = 0))
cat("fitted at k = 1..4, z = 0: Estimate", signif(f[, 1], 6),
    "; Est.Error", signif(f[, 2], 6), "\n")
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(brm(y ~ mo(x), data = dg, chains = 2, iter = 2000, seed = 3,
            refresh = 0))
f <- fitted(fg, newdata = data.frame(x = 1:5))
cat("gap fit, fitted at x = 1..5: Estimate", signif(f[, 1], 6),
    "; Est.Error", signif(f[, 2], 6), "\n")
