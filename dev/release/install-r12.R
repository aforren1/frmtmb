# Roxygenise then install all eight packages into the release library,
# core first: an extension's roxygenise loads the core it will run
# against, so core must be installed before any extension is documented.
lib <- "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb"
pkgs <- c(root, file.path(root, "extensions", paste0("frmtmb.",
  c("eam", "sample", "coupling", "spline", "learn", "latent", "ode"))))
r <- file.path(R.home("bin"), "R")
for (p in pkgs) {
  cat("== roxygenise", basename(p), "\n")
  callr_ok <- tryCatch({
    system2(file.path(R.home("bin"), "Rscript"),
            c("-e", shQuote(sprintf(
              ".libPaths(c('%s','C:/Users/adf44/AppData/Local/R/win-library/4.6')); roxygen2::roxygenise('%s')",
              lib, p))), stdout = TRUE, stderr = TRUE)
  }, error = function(e) conditionMessage(e))
  cat(tail(callr_ok, 3), sep = "\n")
  cat("== install", basename(p), "\n")
  out <- system2(r, c("CMD", "INSTALL", paste0("--library=", lib), shQuote(p)),
                 stdout = TRUE, stderr = TRUE)
  cat(tail(out, 2), sep = "\n")
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("install failed: ", basename(p))
}
for (p in pkgs) {
  nm <- basename(p)
  cat(nm, as.character(packageVersion(nm, lib.loc = lib)), "\n")
}
