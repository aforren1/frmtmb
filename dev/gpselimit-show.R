# Lane gpselimit: what a user reads on one row of the sweep: summary(),
# confint(), the lost reasons and the boundary kinds the check kept.
# Usage: Rscript gpselimit-show.R <lane|base> <row>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
gs_libs(a[1])
row <- as.integer(a[2])
s <- gs_sweep_row(row)
fit <- withCallingHandlers(
  frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian()),
  warning = function(w) {
    cat("WARNING", class(w)[1], ":", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("MESSAGE", class(m)[1], ":", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  })
sdr <- frmtmb:::sdr_of(fit)
cat("se_lost:\n"); print(sdr$se_lost)
cat("boundary kind:\n"); print(fit$cache$se_boundary_kind)
cat("limit info:\n"); str(fit$cache$se_limit)
print(summary(fit))
print(confint(fit))
print(fit$opt$par)
H <- fit$cache$hessian_fixed$H
if (!is.null(H)) {
  cat("outer Hessian:\n"); print(H)
  D <- sqrt(abs(diag(H)))
  print(eigen(H / outer(D, D), symmetric = TRUE))
}
