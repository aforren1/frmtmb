# Lane surface: brms 2.23.0's update.brmsfit, read for item 3.
source("dev/surface-env.R")
surface_env("base")
ns <- asNamespace("brms")
for (f in c("update.brmsfit", "check_prior", ".validate_prior", "validate_prior")) {
  cat("\n####", f, "\n")
  g <- tryCatch(get(f, envir = ns), error = function(e) NULL)
  if (is.null(g)) cat("(not found)\n") else { print(args(g)); print(body(g)) }
}
