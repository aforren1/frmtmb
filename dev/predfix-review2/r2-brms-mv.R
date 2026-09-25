# brms 2.23.0 side of minor 5, a second design: a gaussian, a
# categorical and a second categorical sharing category labels, and a
# cumulative. Seed 606. Compiled with cmdstanr (rstan's stanc.js needs V8).
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
set.seed(606)
n <- 150
d <- data.frame(x = rnorm(n))
d$y1 <- 1 + 0.5 * d$x + rnorm(n)
d$c1 <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = c("a", "b", "c"))
d$c2 <- cut(-d$x + rnorm(n), c(-Inf, 0, Inf), labels = c("a", "b"))
d$o1 <- as.integer(cut(0.8 * d$x + rnorm(n), c(-Inf, -0.7, 0.2, 1, Inf)))
saveRDS(d, "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-brms-mv-data.rds")
b <- brm(mvbf(bf(o1 ~ x, family = cumulative()),
              bf(y1 ~ x, family = gaussian()),
              bf(c1 ~ x, family = categorical()),
              bf(c2 ~ x, family = categorical())), data = d,
         chains = 1, iter = 1000, refresh = 0, seed = 1,
         backend = "cmdstanr")
show <- function(lab, f) {
  cat(lab, ": dim", dim(f), "\n   layers:",
      paste0("'", dimnames(f)[[3]], "'", collapse = " "), "\n")
}
f <- fitted(b)
show("fitted(b)", f)
show("fitted(b, resp = c('c2', 'y1'))", fitted(b, resp = c("c2", "y1")))
show("fitted(b, resp = c('y1', 'o1'))", fitted(b, resp = c("y1", "o1")))
show("fitted(b, resp = 'o1')", fitted(b, resp = "o1"))
nd <- data.frame(x = c(-1, 0, 2))
fn <- fitted(b, newdata = nd)
show("fitted(b, newdata = 3 rows)", fn)
cat("Estimate at newdata:\n")
print(round(fn[, "Estimate", ], 4))
saveRDS(list(f = f, fn = fn), "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-brms-mv.rds")
