# Lane wt-conditions: the evidence behind the class-vector decision.
#   Rscript dev/conditions-design-evidence.R > dev/conditions-log/design.txt
# Reads brms 2.23.0's sources from dev/brms-suite/brms (the CRAN
# tarball, see dev/brms-suite-audit.md) and the installed libraries.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- "dev/brms-suite/brms"
grep_count <- function(dir, pat) {
  fs <- list.files(file.path(src, dir), "[.][Rr]$", full.names = TRUE,
                   recursive = TRUE)
  hits <- unlist(lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    i <- grep(pat, x)
    if (length(i)) paste0(sub(paste0(src, "/"), "", f), ":", i, ": ",
                          trimws(x[i]))
  }))
  cat(sprintf("\n-- %s in brms %s/: %d lines\n", pat, dir, length(hits)))
  if (length(hits) <= 12L) cat(paste0("   ", hits, "\n"), sep = "")
}
cat("== what brms's code and suite rely on ==\n")
for (d in c("R", "tests")) {
  for (p in c("brms_error", "rlang_error", "catch_cnd", "rlang::abort",
              "rlang::warn", "rlang::inform", "cnd_", "trace_back",
              "[^.[:alnum:]_]stop2[(]", "[^.[:alnum:]_]warning2[(]",
              "[^.[:alnum:]_]stopifnot[(]", "[^.[:alnum:]_]stop[(]",
              "[^.[:alnum:]_]warning[(]", "[^.[:alnum:]_]message[(]")) {
    grep_count(d, p)
  }
}
cat("\n== brms's helpers as installed ==\n")
print(brms:::stop2)
print(brms:::warning2)
cat("message2 exists:", exists("message2", asNamespace("brms")), "\n")

cat("\n== the class brms raises ==\n")
e <- tryCatch(brms:::stop2("x"), error = identity)
print(class(e))

cat("\n== is rlang already in frmtmb's dependency closure? ==\n")
db <- installed.packages()
imp <- c("generics", "graphics", "grDevices", "Matrix", "methods", "mgcv",
         "nlme", "reformulas", "RTMB", "RTMBdist", "stats", "TMB", "utils")
clo <- tools::package_dependencies(imp, db = db, recursive = TRUE,
                                   which = c("Depends", "Imports",
                                             "LinkingTo"))
all_deps <- unique(c(imp, unlist(clo)))
cat("closure size:", length(all_deps), "; rlang in closure:",
    "rlang" %in% all_deps, "\n")
cat("closure:", sort(all_deps), "\n")
p <- find.package("rlang")
sz <- sum(file.size(list.files(p, recursive = TRUE, full.names = TRUE)),
          na.rm = TRUE)
cat(sprintf("rlang %s installed size: %.1f MB\n",
            format(packageVersion("rlang")), sz / 2^20))
