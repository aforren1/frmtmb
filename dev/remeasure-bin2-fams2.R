# Lane remeasure, bin 2P G3/G4: which reference densities brms 2.23.0
# exports, which suggested packages the blocks skip on, and the dpars of
# the families frmtmb reaches only through brmsfamily().
#
#   Rscript dev/remeasure-bin2-fams2.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ex <- getNamespaceExports("brms")
for (f in c("dstudent_t", "dshifted_lnorm", "dexgaussian", "dskew_normal",
            "dbeta_binomial", "dfrechet", "dinv_gaussian",
            "dgen_extreme_value", "rasym_laplace", "dcom_poisson",
            "dasym_laplace", "dvon_mises", "pvon_mises", "inv_logit_scaled")) {
  cat(sprintf("%-20s exported by brms: %s\n", f, f %in% ex))
}
for (p in c("extraDistr", "betareg", "RWiener")) {
  cat(p, "installed:", requireNamespace(p, quietly = TRUE), "\n")
}
for (f in c("gaussian", "poisson", "binomial", "gamma", "inverse.gaussian",
            "beta", "com_poisson", "negbinomial2", "frechet",
            "discrete_weibull", "gen_extreme_value")) {
  fam <- tryCatch(brmsfamily(f), error = function(e) e)
  if (inherits(fam, "error")) {
    cat(sprintf("brmsfamily('%s'): ERROR %s\n", f,
                substr(conditionMessage(fam), 1, 60)))
  } else {
    cat(sprintf("brmsfamily('%s'): family=%s dpars=%s type=%s lcdf=%s\n", f,
                fam$family, paste(fam$dpars, collapse = ","), fam$type,
                is.function(fam$lcdf)))
  }
}
