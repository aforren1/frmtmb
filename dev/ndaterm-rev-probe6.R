# Reviewer of lane ndaterm: zero_inflated_binomial at zi near 0. Seed 7.
args <- commandArgs(TRUE); arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE),
                w = runif(n, 0.5, 2), g = gl(10, 20))
d$y <- rbinom(n, d$N, plogis(0.3 + 0.5 * d$x))
fz <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(), data = d)
print(summary(fz)$fixed)
set.seed(1)
p <- predict(fz, ndraws = 400)
cat("in-sample predict mean / fitted mean:", mean(p[, "Estimate"]),
    mean(fitted(fz)[, "Estimate"]), "\n")
