# Parse every workflow this lane changed, and print what the gate steps
# and the matrix read, so the YAML is checked as GitHub would load it
# (keys and values, not only syntax).
#   Rscript dev/smallfix-yaml-parse.R
# Log: dev/smallfix-out/yaml-parse.txt
.libPaths(c(.libPaths(), "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
`%||%` <- function(a, b) if (is.null(a)) b else a
cat("yaml", format(packageVersion("yaml")), "\n")
wf <- c(".github/workflows/R-CMD-check.yaml",
        ".github/workflows/test-coverage.yaml",
        ".github/workflows/absent-helper-identical.yaml",
        ".github/workflows/pkgdown.yaml")
for (f in wf) {
  y <- yaml::read_yaml(f)
  cat("==", f, ": parsed, jobs:", paste(names(y$jobs), collapse = ", "), "\n")
  for (j in names(y$jobs)) {
    jb <- y$jobs[[j]]
    for (cfg in jb$strategy$matrix$config) {
      cat("   matrix:", cfg$os, cfg$r, "absent =", cfg$absent %||% "-", "\n")
    }
    for (st in jb$steps) {
      if (any(grepl("skip-gate|byte-identical|saved-fit-check", c(st$name, st$run)))) {
        cat("   step [", st$name, "] shell:", st$shell %||% "-", "\n")
        cat("     run:", gsub("\n", " | ", trimws(st$run)), "\n")
        if (!is.null(st$env)) cat("     env:", paste(names(st$env), unlist(st$env),
                                                  sep = "=", collapse = " "), "\n")
      }
    }
    sx <- Filter(function(s) identical(s$uses,
                                       "r-lib/actions/setup-r-dependencies@v2"),
                 jb$steps)
    for (s in sx) cat("   setup-r-dependencies extra-packages:",
                      s$with$`extra-packages`, "\n")
  }
}
