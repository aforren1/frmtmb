# Review: fitted(newdata =) of y ~ exp(x) + z + sar(W) (lag and error)
# against (I - rho W)^-1 X_new b computed here; fcor(V) with exp(x).
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spatial)})
k <- 6
g <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A); n <- nrow(W)
set.seed(20261024)
d <- data.frame(x = runif(n, 0, 1.5), z = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * exp(d$x) - 0.3 * d$z +
                          rnorm(n, 0, 0.5)))
nd <- transform(d, x = x + 0.2)
for (ty in c("lag", "error")) {
  fit <- frm(y ~ exp(x) + z + sar(W, type = ty), data = d,
             data2 = list(W = W))
  b <- fixef(fit)[, "Estimate"]
  rho <- unname(fit$frame$autocor[[1]]$aux$lo + (fit$frame$autocor[[1]]$aux$hi -
    fit$frame$autocor[[1]]$aux$lo) * plogis(fit$estimates$thetaac))
  rho2 <- confint_varcorr(fit)
  mu <- as.vector(cbind(1, exp(nd$x), nd$z) %*% b[c("Intercept", "expx", "z")])
  want <- if (ty == "lag") as.vector(solve(diag(n) - rho * W, mu)) else mu
  got <- unname(fitted(fit, newdata = nd)[, "Estimate"])
  cat(ty, ": rownames(fixef) =", rownames(fixef(fit)), "; max rel diff fitted(newdata) =",
      max(abs(got - want)) / max(abs(want)), "\n")
  pr <- unname(predict(fit, newdata = nd, summary = TRUE)[, "Estimate"])
  cat(ty, ": predict(newdata) mean vs fitted:", max(abs(pr - got)) / max(abs(got)), "\n")
}
V <- 0.5^abs(outer(1:n, 1:n, "-"))
fit <- frm(y ~ exp(x) + fcor(V), data = d, data2 = list(V = V))
b <- fixef(fit)[, "Estimate"]
got <- unname(fitted(fit, newdata = nd)[, "Estimate"])
want <- as.vector(cbind(1, exp(nd$x)) %*% b[c("Intercept", "expx")])
cat("fcor: rownames(fixef) =", rownames(fixef(fit)), "; max rel diff =",
    max(abs(got - want)) / max(abs(want)), "\n")
cat("conditional_effects names:", names(conditional_effects(fit)), "\n")
