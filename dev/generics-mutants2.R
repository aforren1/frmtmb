# Four more scale mutants, so that every block of
# tests/testthat/test-scale-contract.R is seen failing against a build
# that reports the wrong scale.  Same shape as dev/generics-mutants.R.
#
#   M4  residuals() default is formed on the LINK scale
#   M5  simulate() returns draws on the LINK scale
#   M6  fitted() returns the median while predict(type = "response")
#       keeps the mean, so the two drift apart
#   M7  VarCorr() reports the covariance on the RESPONSE scale
av <- commandArgs(trailingOnly = TRUE)
which <- av[1]
src <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
out <- file.path("C:/Users/adf44/source/r", paste0("generics-mut-", which))
lib <- paste0(out, "-lib")
unlink(out, recursive = TRUE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (d in c("R", "man", "inst", "src", "data")) {
  if (dir.exists(file.path(src, d))) {
    file.copy(file.path(src, d), out, recursive = TRUE)
  }
}
for (f in c("DESCRIPTION", "NAMESPACE")) {
  file.copy(file.path(src, f), file.path(out, f))
}

sub1 <- function(path, old, new) {
  f <- file.path(out, path)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  n <- length(gregexpr(old, txt, fixed = TRUE)[[1]])
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (n != 1L) stop(n, " matches in ", path)
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], f, useBytes = TRUE)
  cat("mutated", path, "\n")
}

if (which == "M4") {
  sub1("R/predict.R",
    "  mu <- response_mean(fam, dp, av)\n  r <- yv - mu\n",
    "  mu <- response_mean(fam, dp, av)\n  r <- log(yv) - log(mu)\n")
} else if (which == "M5") {
  sub1("R/predict.R",
    "  out <- lapply(out, function(v) sim_restore_type(object, rspec, v))",
    "  out <- lapply(out, function(v) log(sim_restore_type(object, rspec, v)))")
} else if (which == "M6") {
  sub1("R/predict.R",
    paste0("fitted.frmtmb_fit <- function(object, ...) {\n",
           "  rspec <- single_response(object, \"fitted()\")"),
    paste0("fitted.frmtmb_fit <- function(object, ...) {\n",
           "  rspec <- single_response(object, \"fitted()\")\n",
           "  if (identical(rspec$family[[\"family\"]], \"lognormal\")) {\n",
           "    return(exp(predict(object, type = \"link\")))\n",
           "  }"))
} else if (which == "M7") {
  sub1("R/methods-fit.R",
    paste0("  names(out) <- vapply(x$frame[[\"re_blocks\"]], `[[`, \"\",",
           " \"term_label\")\n  structure(out, class = \"VarCorr_frmtmb\")"),
    paste0("  sc <- mean(stats::fitted(x))^2\n",
           "  out <- lapply(out, function(V) V * sc)\n",
           "  names(out) <- vapply(x$frame[[\"re_blocks\"]], `[[`, \"\",",
           " \"term_label\")\n  structure(out, class = \"VarCorr_frmtmb\")"))
} else {
  stop("unknown mutant")
}

o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(lib)),
               "--no-multiarch", shQuote(out)),
             stdout = TRUE, stderr = TRUE)
cat(tail(o, 3), sep = "\n")
cat("LIB", lib, "\n")
