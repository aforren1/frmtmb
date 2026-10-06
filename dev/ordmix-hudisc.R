# Is disc hidden on hurdle_cumulative() under thres(gr = ) when held at
# 1, and shown and fitted when modeled? Seed 20261020.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261020)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- ifelse(runif(n) < 0.2, 0L, 1L + (u > -1) + (u > 0.3) + (u > 1.5))
f0 <- frm(bf(y | thres(gr = g) ~ x), family = hurdle_cumulative(), data = d)
f1 <- frm(bf(y | thres(gr = g) ~ x, disc ~ 0 + z),
          family = hurdle_cumulative(), data = d)
cat("held at 1: links line:", grep("Links", capture.output(print(f0)),
                                    value = TRUE), "\n")
cat("held at 1: fixef(flatten) names:", names(fixef(f0, flatten = TRUE)),
    "\n")
cat("held at 1: variables:", variables(f0), "\n")
cat("modeled: links line:", grep("Links", capture.output(print(f1)),
                                  value = TRUE), "\n")
cat("modeled: disc_z", fixef(f1)["disc_z", ], "\n")
cat("logLik held / modeled:", as.numeric(logLik(f0)), as.numeric(logLik(f1)),
    "\n")
