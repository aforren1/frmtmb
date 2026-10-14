.libPaths(c("C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(31); n <- 200
d <- data.frame(x = rnorm(n)); P <- cbind(1, exp(0.3 + 0.5 * d$x), 1); P <- P / rowSums(P)
d$N <- 10L; d$Ym <- t(vapply(1:n, function(i) rmultinom(1, 10, P[i, ])[, 1], numeric(3)))
f <- frm(bf(Ym | trials(N) ~ x), family = multinomial(K = 3), data = d)
cat(tryCatch({dharma_residuals(f, nsim = 20, seed = 1); "OK"}, error = function(e) conditionMessage(e)), "\n")
