# Review of lane hygiene, item 131: the vignette says a constant weight
# on every row of a model with no random effect "leaves every estimate
# unchanged, sigma included". Is that so under REML = TRUE and with a
# prior on a coefficient? Seed 31, the data of test-weights-constant.R.
#   Rscript dev/reviews/hygiene-review-wreml.R
.libPaths(c("C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(31)
n <- 80L
d <- data.frame(x = rnorm(n), z = runif(n))
d$y <- 1 + 0.6 * d$x - 0.4 * d$z + rnorm(n, 0, 0.7)
d$w2 <- 2
sig <- function(f) exp(f$opt$par[grep("sigma|betad", names(f$opt$par))])
show <- function(label, a, b) {
  pa <- a$opt$par
  pb <- b$opt$par
  cat(sprintf("%-28s max |par diff| %.3g  (%s)\n", label,
              max(abs(pa - pb)),
              paste(sprintf("%s %.5f/%.5f", names(pa), pa, pb),
                    collapse = ", ")))
}
f1 <- frm(y ~ x + z, data = d, REML = TRUE)
f2 <- frm(y | weights(w2) ~ x + z, data = d, REML = TRUE)
show("REML, w = 2 vs none", f1, f2)
cat("REML sigma ratio (w = 2 / none):",
    format(exp(f2$opt$par[["betad"]] - f1$opt$par[["betad"]]), digits = 6),
    "; RSS/(n-p) vs 2 RSS/(2n-p) predicts",
    format(sqrt((2 * (n - 3)) / (2 * n - 3)), digits = 6), "\n")
pr <- set_prior("normal(0, 0.1)", class = "b", coef = "x")
g1 <- frm(y ~ x + z, data = d, prior = pr)
g2 <- frm(y | weights(w2) ~ x + z, data = d, prior = pr)
show("prior on b_x, w = 2 vs none", g1, g2)
h1 <- frm(y ~ x + z, data = d)
h2 <- frm(y | weights(w2) ~ x + z, data = d)
show("ML, w = 2 vs none", h1, h2)
