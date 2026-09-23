source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
cat("bayesplot", format(packageVersion("bayesplot")), "\n")
ty <- as.character(bayesplot::available_ppc(""))
cat(length(ty), "types\n")
ns <- asNamespace("bayesplot")
for (t in ty) {
  f <- names(formals(get(t, ns)))
  cat(sprintf("%-32s %s\n", t, paste(f, collapse = ",")))
}
