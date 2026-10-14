# Reviewer round 2: brms's R-side com_poisson helpers, as installed
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
print(body(brms:::log_lik_com_poisson))
print(body(brms:::posterior_predict_com_poisson))
k <- 0:200
p <- exp(k * log(3) - 2 * lfactorial(k)); p <- p / sum(p)
cat("Stan side mu 3 nu 2 mean", sprintf("%.6f", sum(k * p)), "\n")
