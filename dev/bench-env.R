## environment probe for the optimParallel benchmark
cat("R", R.version.string, "\n")
cat("libs:\n"); print(.libPaths())
for (p in c("frmtmb", "optimParallel", "lme4", "pkgload", "RTMB", "TMB",
            "Matrix", "parallel")) {
  v <- tryCatch(as.character(utils::packageVersion(p)),
                error = function(e) NA_character_)
  cat(sprintf("  %-14s %s\n", p, v))
}
cat("cores (detectCores):", parallel::detectCores(), "\n")
cat("logical/physical:", parallel::detectCores(logical = TRUE), "/",
    parallel::detectCores(logical = FALSE), "\n")
