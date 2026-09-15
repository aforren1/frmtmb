# Item 2.5d's evidence standard: each scale assertion must be SEEN
# FAILING against a build that reports the wrong scale.  Mutating the
# test's own expectation would prove nothing, so this mutates the
# PACKAGE and runs the unmodified test file against it.
#
#   Rscript dev/generics-mutants.R <which>
#
#   M1  lognormal fitted() returns the MEDIAN exp(mu) instead of the
#       mean exp(mu + sigma^2 / 2)
#   M2  predict() defaults to the RESPONSE scale, as brms does
#   M3  sigma() returns the LINK scale (log sigma), not the response
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

if (which == "M1") {
  sub1("R/families.R",
    paste0("      mean_fn = function(dpars, aterms) {\n",
           "        exp(dpars[[\"mu\"]] + dpars[[\"sigma\"]]^2 / 2)\n",
           "      },"),
    paste0("      mean_fn = function(dpars, aterms) {\n",
           "        exp(dpars[[\"mu\"]])\n",
           "      },"))
} else if (which == "M2") {
  sub1("R/predict.R",
    paste0("                               type = c(\"link\", \"response\",\n",
           "                                        \"conditional\", \"zprob\", \"zlink\",\n",
           "                                        \"disp\"),"),
    paste0("                               type = c(\"response\", \"link\",\n",
           "                                        \"conditional\", \"zprob\", \"zlink\",\n",
           "                                        \"disp\"),"))
} else if (which == "M3") {
  sub1("R/sugar.R",
    "      return(lp[[\"link\"]]$linkinv(object$estimates[[lp[[\"par\"]]]][lp[[\"idx\"]]]))",
    "      return(object$estimates[[lp[[\"par\"]]]][lp[[\"idx\"]]])")
} else {
  stop("unknown mutant")
}

o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(lib)),
               "--no-multiarch", shQuote(out)),
             stdout = TRUE, stderr = TRUE)
cat(tail(o, 3), sep = "\n")
cat("LIB", lib, "\n")
