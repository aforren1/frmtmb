# Lane optstart: the ported brms fixture count ~ Trt * Age + (1 +
# mmc(Age, volume) | mm(patient, visit)) on brmsfit_example3's data, where
# the Hessian scan found the lane's outer Hessian not positive definite
# and 0.70.0's positive definite. Objective, code, units, warnings, the
# SEs vcov() gives and a Richardson Jacobian's.
#   Rscript dev/optstart-mm.R <lib or "">
.libPaths(c(if (nzchar(commandArgs(TRUE)[1])) commandArgs(TRUE)[1],
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
# the data of brms:::brmsfit_example3, as helper-brms-suite.R reads them
d <- as.data.frame(get("brmsfit_example3", envir = asNamespace("brms"))$data)
attr(d, "terms") <- NULL
w <- character()
f <- withCallingHandlers(suppressMessages(frm(
  bf(count ~ Trt * Age + (1 + mmc(Age, volume) | mm(patient, visit))),
  family = gaussian(), data = d)), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
cat(sprintf("objective %.10f code %d | units %s\n", f$opt$objective,
            f$opt$convergence, paste(signif(f$par_units %||% 1, 3),
                                     collapse = " ")))
se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
w2 <- character()
withCallingHandlers(invisible(vcov(f)), warning = function(x) {
  w2 <<- c(w2, conditionMessage(x)); invokeRestart("muffleWarning")
})
H <- numDeriv::jacobian(function(x) as.numeric(f$obj$gr(x)), f$opt$par)
ref <- sqrt(diag(solve((H + t(H)) / 2)))
print(signif(cbind(vcov = se, richardson = ref, rel = se / ref - 1), 4))
cat("warnings at fit:", substr(paste(w, collapse = " | "), 1, 200), "\n")
cat("warnings at vcov:", substr(paste(w2, collapse = " | "), 1, 200), "\n")
