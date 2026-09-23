# Punch round 2: the numbers behind each fix, on the round-2 build.
#   Rscript dev/shapes-p2-measure.R > dev/shapes-log/p2-measure.txt
.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260921)
n <- 150
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
u <- rnorm(15, 0, 0.8)
d$y <- rnorm(n, 1 + 0.5 * d$x + u[d$g], 1)
e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2))
d$y1 <- 1 + 0.5 * d$x + e[, 1]
d$y2 <- -0.3 * d$x + e[, 2]
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))

cat("== A. unseen levels, propagate_error = FALSE, 4000 draws\n")
fm <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
nd <- data.frame(x = 0, g = factor(c("new1", "new2", "new3", "new4", "new4")))
set.seed(1)
dr <- predict(fm, newdata = nd, allow_new_levels = TRUE,
              propagate_error = FALSE, ndraws = 4000, summary = FALSE)
tau2 <- VarCorr(fm)$g$sd[1L, "Estimate"]^2
s2 <- sigma(fm)^2
print(round(cor(dr), 4))
cat("tau^2/(tau^2+sigma^2) =", round(tau2 / (tau2 + s2), 4), "\n\n")

cat("== B. poisson row alone and beside an overflowing row, 400 draws\n")
fp <- frm(bf(cnt ~ x) + poisson(), data = d)
b <- fixef_by_dpar(fp)$mu
x_bad <- (709.5 - b[["(Intercept)"]]) / b[["x"]]
set.seed(7)
a1 <- predict(fp, newdata = data.frame(x = 4), ndraws = 400)
set.seed(7)
a2 <- withCallingHandlers(
  predict(fp, newdata = data.frame(x = c(4, x_bad)), ndraws = 400),
  warning = function(w) {
    cat("warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
print(rbind(alone = a1[1, ], beside = a2[1, ]))
cat("identical:", identical(a1[1, ], a2[1, ]), "\n\n")

cat("== C. rescor draws, 3000 draws\n")
fr <- frm(bf(mvbind(y1, y2) ~ x) + gaussian() + set_rescor(TRUE), data = d)
for (pu in c(FALSE, TRUE)) {
  set.seed(3)
  dd <- predict(fr, summary = FALSE, ndraws = 3000, propagate_error = pu)
  w <- vapply(seq_len(dim(dd)[2L]), function(i) cor(dd[, i, 1], dd[, i, 2]), 0)
  cat("propagate_error =", pu, " mean within-row draw correlation",
      round(mean(w), 4), "\n")
}
cat("estimated rescor", round(rescor_matrix(fr)[1, 2], 4), "\n")
