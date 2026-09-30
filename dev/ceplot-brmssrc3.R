# Third dump of brms 2.23.0 internals for lane ceplot: small helpers the
# plot methods call.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("brms")
for (nm in c("use_alias", "limit_chars", "replace_args<-", "is.theme",
             "make_point_frame", "conditional_smooths.brmsfit")) {
  cat("#####", nm, "\n")
  f <- get(nm, envir = ns)
  if (is.function(f)) print(f) else print(f)
}
