# lane tmbstan, probe 01: the installed tmbstan/StanHeaders state.
# Reproduces the measurement quoted in the lane brief. Read-only.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

hpp <- system.file("model.hpp", package = "tmbstan")
cat("tmbstan libpath:", find.package("tmbstan"), "\n")
cat("tmbstan version:",
    as.character(utils::packageVersion("tmbstan")), "\n")
cat("hpp found:", nzchar(hpp), "\n")
if (nzchar(hpp)) {
  ln <- readLines(hpp, warn = FALSE)
  cat("hpp lines:", length(ln), "\n")
  marker <- "std_normal_lpdf<propto__>(y)"
  hit <- grep(marker, ln, fixed = TRUE)
  cat("broken marker present:", length(hit) > 0, "\n")
  cat("marker line numbers:", paste(hit, collapse = ","), "\n")
  # how many log_prob_impl overloads the generator emitted, which is
  # the upstream fact the guard is a proxy for
  cat("log_prob_impl occurrences:",
      length(grep("log_prob_impl", ln, fixed = TRUE)), "\n")
  cat("std_normal_lpdf occurrences (any form):",
      length(grep("std_normal_lpdf", ln, fixed = TRUE)), "\n")
  patched <- grep("objective_function", ln, fixed = TRUE)
  cat("objective_function occurrences:", length(patched), "\n")
}
cat("StanHeaders:",
    as.character(utils::packageVersion("StanHeaders")), "\n")
cat("StanHeaders libpath:", find.package("StanHeaders"), "\n")
cat("rstan:", as.character(utils::packageVersion("rstan")), "\n")
# build-time evidence: when was tmbstan installed, and against what
d <- read.dcf(file.path(find.package("tmbstan"), "DESCRIPTION"))
for (f in intersect(c("Built", "Packaged", "Date/Publication"),
                    colnames(d))) {
  cat(f, ": ", d[, f], "\n", sep = "")
}
