# Reviewer: a new group's variance on y ~ mo(x) * f2 + (mo(x) | g):
# se^2 with allow_new_levels minus se^2 at re_formula = NA should be
# (1, m(x)) S (1, m(x))', m the term's transform at the row's category.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 90
d <- data.frame(x = sample(0:3, n, TRUE), f2 = factor(rep(c("a", "b"), 45)),
                g = factor(rep(1:9, 10)))
d$y <- 1 + 0.5 * c(0, 0.5, 0.7, 1.4)[d$x + 1] * (1 + rnorm(9, 0, 0.4)[d$g]) +
  rnorm(9, 0, 0.4)[d$g] + rnorm(n, 0, 0.5)
fr <- suppressWarnings(frm(bf(y ~ mo(x) * f2 + (mo(x) | g)),
                           family = gaussian(), data = d))
nd <- data.frame(x = 0:3, f2 = factor("a", levels = c("a", "b")),
                 g = factor("new"))
a <- frm_linpred(fr, newdata = nd, se.fit = TRUE, allow_new_levels = TRUE)
b <- frm_linpred(fr, newdata = nd, se.fit = TRUE, re_formula = NA)
lp <- fr$frame$linpreds[[1]]
mi <- lp$mo[[1]]
m <- frmtmb:::mo_col_values(fr, mi, 0:3)
vc <- VarCorr(fr)$g
S <- vc$cov[, 1, ]
hand <- vapply(1:4, function(i) {
  v <- c(1, m[i])
  as.numeric(t(v) %*% S %*% v)
}, 0)
cat("transform m(x):", format(m, digits = 5), "\n")
cat("se^2 new - se^2 NA:", format(a$se.fit^2 - b$se.fit^2, digits = 6), "\n")
cat("hand (1, m) S (1, m)':", format(hand, digits = 6), "\n")
cat("eta new - eta NA:", format(a$fit - b$fit, digits = 4), "\n")
