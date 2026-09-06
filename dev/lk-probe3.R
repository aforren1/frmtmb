for (f in c("softit","inv_softit","log_expm1","log1p_exp","logit","inv_logit","cloglog","inv_cloglog","expp1","logm1","incl_dpars")) {
  cat("=====", f, "=====\n")
  ok <- tryCatch({print(get(f, envir = asNamespace("brms"))); TRUE}, error = function(e) FALSE)
  if (!ok) cat("(absent)\n")
}
