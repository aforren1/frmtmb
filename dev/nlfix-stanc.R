# Which of brms 2.23.0's names_stan_functions() are functions of the Stan
# language a brms body can call, as rstan 2.32 parses them? stanc only, no
# compile. Signature tried: real arguments, one or two, then vector.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(rstan))
nm <- brms:::names_stan_functions()
try_sig <- function(f, call_txt, type = "real") {
  code <- sprintf(paste0("data { real a; real b; vector[3] v; }\n",
                         "generated quantities { %s r = %s; }\n"), type,
                  call_txt)
  ok <- tryCatch({
    suppressMessages(capture.output(rstan::stanc(model_code = code)))
    TRUE
  }, error = function(e) FALSE)
  ok
}
res <- lapply(nm, function(f) {
  a <- c(
    real1 = try_sig(f, sprintf("%s(a)", f)),
    real2 = try_sig(f, sprintf("%s(a, b)", f)),
    vec1 = try_sig(f, sprintf("%s(v)", f), "vector[3]"),
    vec_to_real = try_sig(f, sprintf("%s(v)", f)))
  data.frame(name = f, t(a))
})
tab <- do.call(rbind, res)
print(tab, row.names = FALSE)
cat("in Stan:", paste(tab$name[rowSums(tab[, -1]) > 0], collapse = ", "), "\n")
cat("not in Stan:", paste(tab$name[rowSums(tab[, -1]) == 0], collapse = ", "),
    "\n")
