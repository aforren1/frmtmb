# Ordinal mixture: predict() shares against fitted(), with and without
# propagate_error, to tell a simulator defect from a draw-law effect.
# Usage: Rscript dev/ppdraw-probe2.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(7)
d <- data.frame(x = stats::rnorm(400))
k <- stats::rbinom(400, 1, 0.3)
d$yo <- cut(d$x + ifelse(k == 1, 2, 0) + stats::rlogis(400),
            c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
fit <- frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
           data = d)
nd <- data.frame(x = c(-1, 0, 1))
fe <- fitted(fit, newdata = nd)[, "Estimate", ]
cat("fitted:\n"); print(round(fe, 4))
B <- 20000L
sh <- function(p) t(apply(p, 2L, function(col) tabulate(col, 4L) / B))
set.seed(13)
p0 <- predict(fit, newdata = nd, summary = FALSE, ndraws = B,
              propagate_error = FALSE)
s0 <- sh(p0)
cat("shares, propagate_error = FALSE:\n"); print(round(s0, 4))
cat("z:\n"); print(round((s0 - fe) / sqrt(fe * (1 - fe) / B), 2))
set.seed(13)
p1 <- suppressMessages(predict(fit, newdata = nd, summary = FALSE,
                               ndraws = B))
s1 <- sh(p1)
cat("shares, propagate_error = TRUE:\n"); print(round(s1, 4))
cat("z:\n"); print(round((s1 - fe) / sqrt(fe * (1 - fe) / B), 2))
# the mean of the category probabilities over the parameter draws, by
# hand, from the same law predict() draws from
ds <- frmtmb:::fit_draw_space(fit)
cat("theta1 estimate", fit$estimates$betad[3], "sd", sqrt(ds$V[3, 3]), "\n")
cat("summary of fit:\n"); print(summary(fit))
