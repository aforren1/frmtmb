# Lane surface: print the brms 2.23.0 functions whose semantics the
# lane matches (read, not run).
#
#   Rscript dev/surface-brms-src.R > dev/surface-out/brms-src.txt
source("dev/surface-env.R")
surface_env("base")
ns <- asNamespace("brms")
for (f in c("add_criterion.brmsfit", "summary.brmsfit", "plot.brmsfit",
            "pp_mixture.brmsfit", "stancode.brmsfit", "standata.brmsfit",
            "fixef.brmsfit", "loo.brmsfit", "compute_loo",
            "compute_loolist")) {
  cat("\n####", f, "\n")
  g <- tryCatch(get(f, envir = ns), error = function(e) NULL)
  if (is.null(g)) {
    cat("(not found)\n")
  } else {
    print(args(g))
    print(body(g))
  }
}
