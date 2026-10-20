# Consolidation item 3: parse .github/workflows/pkgdown.yaml, check every
# Rscript step parses as R, and print the new step; then run the new
# step's R code locally, in this tree, against the library named (all
# seven extensions installed there), as the job would. Adapted from
# dev/testiso-yaml.R.
#   Rscript dev/rel072-yaml.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
y <- yaml::read_yaml(".github/workflows/pkgdown.yaml")
cat("jobs:", names(y$jobs), "\n")
step <- NULL
for (j in names(y$jobs)) {
  for (s in y$jobs[[j]]$steps) {
    if (identical(s$shell, "Rscript {0}")) {
      ok <- tryCatch({ parse(text = s$run); "parses" },
                     error = function(e) conditionMessage(e))
      cat(j, "/", s$name, ":", ok, "\n")
    }
    if (identical(s$name, "Check the extensions' compatibility rows")) {
      step <- s
    }
  }
}
stopifnot(!is.null(step))
cat("new step env:", paste0(names(step$env), "=", unlist(step$env)), "\n")
names_in_job <- vapply(y$jobs[[1]]$steps, function(s) s$name %||% "", "")
cat("step order:", match(c("Install the seven extensions from this checkout",
                            step$name), names_in_job), "\n")
do.call(Sys.setenv, step$env)
cat("---- running the step's code\n")
eval(parse(text = step$run), envir = globalenv())
cat("STEP OK\n")
