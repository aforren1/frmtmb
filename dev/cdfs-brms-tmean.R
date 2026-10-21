# brms 2.23.0's discrete truncated mean against its own likelihood's
# window: posterior_epred_trunc_discrete() at lambda = 3, lb = 2, ub = 6
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
m <- suppressMessages(brms:::posterior_epred_trunc_discrete(
  "pois", list(lambda = matrix(3, 2, 2)), lb = matrix(2, 2, 2),
  ub = matrix(6, 2, 2)))[1, 1]
k_in <- 2:6; w_in <- dpois(k_in, 3)
k_ex <- 3:6; w_ex <- dpois(k_ex, 3)
cat(sprintf("brms epred %.10f  inclusive [2, 6] %.10f  exclusive (2, 6] %.10f\n",
            m, sum(k_in * w_in) / sum(w_in), sum(k_ex * w_ex) / sum(w_ex)))
# the likelihood's window, from brms's Stan code
d <- data.frame(y = c(2, 3, 4, 6))
cat(grep("lccdf|lcdf", strsplit(stancode(y | trunc(lb = 2, ub = 6) ~ 1,
  data = d, family = poisson()), "\n")[[1]], value = TRUE), sep = "\n")
