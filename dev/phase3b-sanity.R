# Quick sanity of the new Wiener distribution function against RWiener,
# during development. Not evidence; the evidence scripts are
# phase3b-cdf-reference.R and phase3b-cdf-validate.R.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE,
                                   export_all = TRUE))
pts <- expand.grid(t = c(0.05, 0.3, 1, 3), v = c(-2, 0, 1.5), a = c(0.8, 2),
                   w = c(0.3, 0.6))
r <- ddm_rt_lcdf2(pts$t, pts$v, pts$a, pts$w)
ref <- mapply(function(t, v, a, w) {
  RWiener::pwiener(t + 1e-9, a, 1e-9, w, v, resp = "both")
}, pts$t, pts$v, pts$a, pts$w)
print(cbind(pts, F = exp(r$lF), ref, S = exp(r$lS),
            sum = exp(r$lF) + exp(r$lS)))
cat("max abs vs RWiener:", max(abs(exp(r$lF) - ref)), "\n")
