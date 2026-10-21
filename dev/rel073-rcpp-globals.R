# Where the Rcpp module classes of a Stan model land (brief item 6). A
# gated one-session run lists `.__C__Rcpp_stan_fit4model<hash>_` and
# its neighbors as globals a test file added. This loads no frmtmb
# package and no test: it reads two cached brms programs from
# dev/stan-cache and instantiates each model's module from inside a
# function, as rstan::sampling() does, then lists what the global
# environment gained. If rstan alone puts the names there, the LEAK line
# is rstan's and run-session.R may tag it.
#
#   Rscript dev/rel073-rcpp-globals.R
#   (log: dev/rel073-log/rcpp-globals.txt)
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
cache <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/stan-cache"
g0 <- ls(globalenv(), all.names = TRUE)
cat("frmtmb loaded:", "frmtmb" %in% loadedNamespaces(), "\n")
suppressMessages(requireNamespace("rstan"))
cat("rstan", format(packageVersion("rstan")), "loaded; globals added by",
    "the load:", setdiff(ls(globalenv(), all.names = TRUE), g0), "\n")
inst <- function(sm) {
  # the step of rstan::sampling() that builds the model's C++ module
  mod <- sm@mk_cppmodule(sm)
  invisible(class(mod))
}
n <- 0L
for (f in sort(list.files(cache, "[.]rds$", full.names = TRUE))) {
  sm <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!inherits(sm, "stanmodel")) next
  before <- ls(globalenv(), all.names = TRUE)
  ok <- tryCatch({ inst(sm); TRUE }, error = function(e) {
    cat("  ", basename(f), "error:", conditionMessage(e), "\n")
    FALSE
  })
  if (!ok) next
  added <- setdiff(ls(globalenv(), all.names = TRUE), before)
  cat(basename(f), "model", sm@model_name, "| globals added:",
      if (length(added)) added else "(none)", "\n")
  n <- n + 1L
  if (n == 2L) break
}
cat("search path:", search(), "\n")
cat("END\n")
