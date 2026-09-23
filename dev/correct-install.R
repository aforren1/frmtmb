# Lane wt-correct: roxygenise and install the touched packages into the
# lane library only. Usage: Rscript dev/correct-install.R [core] [sample]
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
LIB <- Sys.getenv("CORRECT_LIB", "C:/Users/adf44/source/r/correct-lib")
ROOT <- "C:/Users/adf44/source/r/frmtmb-wt-correct"
which <- commandArgs(trailingOnly = TRUE)
dirs <- c(core = ROOT, sample = file.path(ROOT, "extensions/frmtmb.sample"))
for (w in which) {
  d <- dirs[[w]]
  cat("== roxygenise", d, "\n")
  roxygen2::roxygenise(d)
  cat("== install", d, "\n")
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
                  paste0("--library=", LIB), shQuote(d)))
  if (st != 0) stop("install failed: ", d)
}
cat("installed:", vapply(c("frmtmb", "frmtmb.sample"), function(p) {
  tryCatch(paste(p, format(packageVersion(p, lib.loc = LIB))),
           error = function(e) paste(p, "absent"))
}, ""), "\n")
