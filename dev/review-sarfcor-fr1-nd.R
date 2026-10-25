# Reviewer, fix round 1, F3: same-row newdata in frmtmb, error form,
# refusals, and a fit that dropped NA rows.
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib", "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
set.seed(4); k <- 6; g <- expand.grid(r = 1:k, c = 1:k); n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j) as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
dg <- 1 / sqrt(rowSums(A)); W <- A * outer(dg, dg)
d <- data.frame(x = rnorm(n)); d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n)))
nd <- d; nd$x <- nd$x + 1
for (type in c("lag", "error")) {
  f <- frm(as.formula(paste0("y ~ x + sar(W, type = '", type, "')")), data = d, data2 = list(W = W))
  rho <- confint_varcorr(f)$estimate
  b <- fixef(f)[, "Estimate"]; mun <- as.vector(cbind(1, nd$x) %*% b)
  ref <- if (type == "lag") solve(diag(n) - rho * W, mun) else mun
  cat(sprintf("%s: fitted(newdata) vs reference %.3g; predict mean %.3g; simulate dim %s\n", type,
              max(abs(fitted(f, newdata = nd)[, "Estimate"] - ref)),
              max(abs(predict(f, newdata = nd)[, "Estimate"] - ref)),
              paste(dim(simulate(f, nsim = 3, newdata = nd)), collapse = "x")))
  r <- tryCatch({fitted(f, newdata = nd[1:5, ]); "answered"}, error = function(e) conditionMessage(e))
  cat("  5 rows:", substr(r, 1, 160), "\n")
}
dn <- d; dn$x[c(3, 9)] <- NA
f <- frm(y ~ x + sar(W), data = dn, data2 = list(W = W))
cat("NA fit nobs", nobs(f), "\n")
for (lab in c("newdata 36 rows (NAs filled)", "newdata 34 rows (kept)")) {
  ndd <- if (grepl("36", lab)) d else d[-c(3, 9), ]
  r <- tryCatch({v <- fitted(f, newdata = ndd); paste("answered", nrow(v), "rows")},
                error = function(e) conditionMessage(e))
  cat(lab, ":", substr(r, 1, 200), "\n")
}
