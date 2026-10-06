# Reviewer of lane fixes, re-check: check_nl_identified() on identified
# models in natural units, where the guard's fixed points (coefficients
# near 0.15 + 0.35 sin(.)) put every row of a body in saturation, and on
# tiny n.
#   Rscript dev/fixes-rev2-nl.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
run <- function(lab, expr) {
  r <- tryCatch({
    fit <- suppressMessages(suppressWarnings(expr))
    fe <- fixef(fit)
    sprintf("FIT conv=%d est %s se %s", fit$opt$convergence,
            paste(sprintf("%.5g", fe[, "Estimate"]), collapse = " "),
            paste(sprintf("%.3g", fe[, "Est.Error"]), collapse = " "))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 170)))
  cat(sprintf("%-50s %s\n", lab, r))
}
# 1. psychometric function, stimulus in natural units (50 to 150)
set.seed(21)
x <- runif(600, 50, 150)
d1 <- data.frame(x = x, y = rbinom(600, 1, pnorm(x, 100, 12)))
run("psychometric pnorm(x, m0, exp(ls)), x 50..150",
    frm(bf(y ~ pnorm(x, m0, exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
        family = bernoulli(link = "identity"), data = d1,
        start = list(beta = c(100, log(12)))))
run("same, glm probit reference",
    {g <- glm(y ~ x, family = binomial("probit"), data = d1)
     cat("    glm probit: m0 =", -coef(g)[1] / coef(g)[2],
         " sd =", 1 / coef(g)[2], "\n")
     frm(bf(y ~ pnorm((x - m0) / exp(ls)), m0 ~ 1, ls ~ 1, nl = TRUE),
         family = bernoulli(link = "identity"), data = d1,
         start = list(beta = c(100, log(12))))})
# 2. logistic growth over calendar years (SSlogis shape)
set.seed(22)
yr <- seq(1900, 2000, length.out = 120)
d2 <- data.frame(yr = yr, y = 50 / (1 + exp((1950 - yr) / 12)) +
                   rnorm(120, 0, 1.5))
run("logistic Asym / (1 + exp((xmid - yr) / scal))",
    frm(bf(y ~ Asym / (1 + exp((xmid - yr) / scal)), Asym ~ 1, xmid ~ 1,
           scal ~ 1, nl = TRUE), data = d2,
        start = list(beta = c(50, 1950, 12))))
n2 <- nls(y ~ SSlogis(yr, Asym, xmid, scal), data = d2)
cat("    nls:", format(coef(n2), digits = 6), "\n")
# 3. exponential decay, time in seconds 1000 to 5000
set.seed(23)
tt <- runif(150, 1000, 5000)
d3 <- data.frame(tt = tt, y = 10 * exp(-5e-4 * tt) + rnorm(150, 0, 0.2))
run("decay y0 * exp(-exp(lk) * tt), tt 1000..5000",
    frm(bf(y ~ y0 * exp(-exp(lk) * tt), y0 ~ 1, lk ~ 1, nl = TRUE),
        data = d3, start = list(beta = c(10, log(5e-4)))))
# 4. a linear sum, identified, with a covariate in large units
set.seed(24)
d4 <- data.frame(x = runif(100, 1e3, 2e3), z = rnorm(100))
d4$y <- 1 + 0.002 * d4$x + 0.5 * d4$z + rnorm(100, 0, 0.3)
run("a + b, a ~ 1 + x (x 1e3..2e3), b ~ 0 + z",
    frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + z, nl = TRUE), data = d4))
# 5. tiny n: 4 rows, 2 coefficients
d5 <- data.frame(x = c(1, 2, 3, 4), y = c(2.1, 3.9, 8.2, 15.8))
run("tiny n = 4, a * exp(b * x)",
    frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d5,
        start = list(beta = c(1, 0.7))))
# 6. an absent factor level after subsetting, a ~ 0 + f
set.seed(26)
d6 <- data.frame(f = factor(sample(c("p", "q", "r"), 90, TRUE)),
                 x = runif(90))
d6$y <- c(p = 1, q = 2, r = 3)[as.character(d6$f)] * exp(0.4 * d6$x) +
  rnorm(90, 0, 0.1)
d6 <- d6[d6$f != "r", ]
run("absent level r, a ~ 0 + f, a * exp(b * x)",
    frm(bf(y ~ a * exp(b * x), a ~ 0 + f, b ~ 1, nl = TRUE), data = d6,
        start = list(beta = c(1, 1, 0))))
