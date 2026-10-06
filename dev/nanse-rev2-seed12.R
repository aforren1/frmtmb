# Reviewer, punch round 1: mo() seed 12, whose conditional_effects()
# printed in round 0 and fails to print now. What is lost, why, and are
# the NaN bands honest?
#   Rscript dev/nanse-rev2-seed12.R
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
grDevices::pdf(NULL)
set.seed(12)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
d <- data.frame(income, ls)
d$age <- rnorm(100, mean = 40, sd = 10)
w <- character()
fit <- withCallingHandlers(frm(ls ~ mo(income) * age, data = d),
  warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
for (x in w) cat("frm warn:", substr(x, 1, 300), "\n")
nm <- ns$outer_par_names(fit)
p <- fit$opt$par
cat("par:", paste(sprintf("%s=%.4g", nm, p), collapse = " "), "\n")
cat("grad:", paste(sprintf("%.2g", fit$obj$gr(p)), collapse = " "), "\n")
H <- fit$obj$he(p)
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("unit-diag eigenvalues:", paste(signif(e$values, 3), collapse = " "),
    "\n")
k <- which.min(e$values)
cat("min direction:", paste(sprintf("%s=%.3f", nm, e$vectors[, k])[
  abs(e$vectors[, k]) > 0.01], collapse = " "), "\n")
lost <- ns$sdr_of(fit)$se_lost
cat("lost:", paste(sprintf("%s(%s)", names(lost), lost), collapse = " "), "\n")
wc <- character()
ce <- withCallingHandlers(conditional_effects(fit, "income:age"),
  warning = function(x) {
    wc <<- c(wc, conditionMessage(x)); invokeRestart("muffleWarning")
  })
cat("ce rows", nrow(ce[[1]]), "finite se__", sum(is.finite(ce[[1]]$se__)),
    "| warnings:", paste(substr(wc, 1, 120), collapse = " || "), "\n")
r <- tryCatch({print(ce); "printed"}, error = function(e) conditionMessage(e))
cat("print:", r, "\n")
# a main-effect display: the lost direction loads 0.022 on zeta1_2,
# which keeps its own SE
for (ef in c("income", "age")) {
  wc <- character()
  ce <- withCallingHandlers(conditional_effects(fit, ef),
    warning = function(x) {
      wc <<- c(wc, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  cat(ef, ": rows", nrow(ce[[1]]), "finite se__", sum(is.finite(ce[[1]]$se__)),
      "| warnings", length(wc), "\n")
}
s <- summary(fit)
print(fixef(fit))
