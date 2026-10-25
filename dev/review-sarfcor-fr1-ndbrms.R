# Reviewer: brms 2.23.0 newdata row counts after a fit that dropped NA rows.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(4); k <- 6; g <- expand.grid(r = 1:k, c = 1:k); n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j) as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
dg <- 1 / sqrt(rowSums(A)); W <- A * outer(dg, dg)
d <- data.frame(x = rnorm(n), y = rnorm(n)); dn <- d; dn$x[c(3, 9)] <- NA
fs <- suppressWarnings(brm(y ~ x + sar(W), data = dn, data2 = list(W = W), empty = TRUE))
for (lab in c("36", "34")) {
  ndd <- if (lab == "36") d else d[-c(3, 9), ]
  r <- tryCatch(paste(dim(standata(fs, newdata = ndd)$Msar), collapse = " x "),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("brms newdata", lab, "rows:", r, "\n")
}
