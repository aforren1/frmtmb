# Item 3: predict(newdata = ) with NA responses under ar(cov = FALSE).
# Seed 31. FORMROBUST_LIB="" for the before arm.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(31)
G <- 30; Tn <- 8
d <- expand.grid(t = 1:Tn, g = factor(1:G))
d$x <- rnorm(nrow(d))
e <- as.vector(apply(matrix(rnorm(G * Tn), Tn, G), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$y <- 1 + 0.5 * d$x + e
fit <- frm(bf(y ~ x + ar(t, g, p = 1)), data = d)
ar1 <- unname(fit$estimates[["thetaac"]][1])
sig <- unname(sigma(fit)[1])
b <- fixef(fit)[, "Estimate"]
cat("ar =", format(ar1, digits = 10), " sigma =", format(sig, digits = 10),
    "\n")
nd <- d[d$g %in% c("1", "2"), ]
nd$y[nd$g == "1" & nd$t >= 5] <- NA   # rows 5..8 of group 1 unknown
r <- tryCatch({
  set.seed(1)
  p <- predict(fit, newdata = nd, ndraws = 4000, propagate_error = FALSE)
  p
}, error = function(e) conditionMessage(e))
if (is.character(r)) {
  cat("predict: ERROR:", r, "\n")
} else {
  mu <- b[["Intercept"]] + b[["x"]] * nd$x
  i5 <- which(nd$g == "1" & nd$t == 5)
  i6 <- which(nd$g == "1" & nd$t == 6)
  # row 5 reads the OBSERVED residual of row 4; row 6 reads a filled one
  m5 <- mu[i5] + ar1 * (nd$y[i5 - 1] - mu[i5 - 1])
  sd6 <- sig * sqrt(1 + ar1^2)
  cat("row 5 Estimate / one-step mean:",
      format(r[i5, "Estimate"] / m5, digits = 6), "\n")
  cat("row 5 Est.Error / sigma:", format(r[i5, "Est.Error"] / sig, digits = 6),
      "\n")
  cat("row 6 Est.Error / (sigma sqrt(1 + ar^2)):",
      format(r[i6, "Est.Error"] / sd6, digits = 6), "\n")
  cat("row 6 Est.Error / sigma:", format(r[i6, "Est.Error"] / sig,
                                         digits = 6), "\n")
  cat("anyNA(Estimate):", anyNA(r[, "Estimate"]), "\n")
  f <- fitted(fit, newdata = nd)
  m6 <- mu[i6] + ar1 * (m5 - mu[i5])
  cat("fitted row 6 / expected-fill mean:", format(f[i6, "Estimate"] / m6,
                                                   digits = 15), "\n")
  nd2 <- nd; nd2$y <- NULL
  set.seed(1)
  p2 <- predict(fit, newdata = nd2, ndraws = 50)
  cat("no response column: anyNA =", anyNA(p2[, "Estimate"]), " dim",
      dim(p2), "\n")
}
