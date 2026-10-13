# Lane optrap: the ported suite's fixture 1 (brms's brmsfit_example1,
# student with mo(), s(), arma and sigma ~ Trt), and the row
# brmsfit-methods:326 (fitted(fit1, dpar = "sigma") all positive).
#   Rscript dev/optrap-fix1.R <lib or "">
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
fit <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit$data)
attr(d, "terms") <- NULL
w <- character()
f <- withCallingHandlers(
  frm(bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
           (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt),
      family = student(), data = d, verbose = 1),
  warning = function(x) {w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")})
cat(sprintf("code %d (%s) obj %.10f evals %s\n", f$opt$convergence,
            f$opt$message, f$opt$objective, f$opt$evals))
cat("units", signif(f$par_units, 3), "\n")
print(round(f$opt$par, 4))
for (x in w) cat("W:", substr(x, 1, 150), "\n")
fi <- fitted(f, dpar = "sigma")
cat("fitted sigma: any NA", anyNA(fi), "all > 0", isTRUE(all(fi > 0)), "\n")
print(head(fi, 3))
lost <- frmtmb:::sdr_of(f)$se_lost
cat("se lost:", paste(names(lost), lost, sep = ":", collapse = ", "), "\n")
print(summary(f)$coefficients$sigma)
