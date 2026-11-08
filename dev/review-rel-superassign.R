# Release review 0.80.0: the superassign reviewer's globals check on the
# release tree (base 3cd2a0c8 from git, not the main checkout's files),
# then every new free variable looked up in the installed release
# namespace or base R.
#   Rscript dev/review-rel-superassign.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
REL <- "C:/Users/adf44/source/r/frmtmb-wt-rel080"
base <- file.path(normalizePath(tempdir(), winslash = "/"), "base")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-superassign/dev/review-superassign-globals.R")
src <- sub('^base_root <- .*$', paste0('base_root <- "', base, '"'), src)
src <- sub('^lane_root <- .*$', paste0('lane_root <- "', REL, '"'), src)
fl <- eval(parse(text = src[grep("^files <- c[(]", src):
                              (grep("^# collect every", src) - 1L)]))
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
for (f in fl) {
  dir.create(dirname(file.path(base, f)), recursive = TRUE, showWarnings = FALSE)
  writeLines(system2("git", c("-C", REL, "-c", "core.autocrlf=false", "show",
                              paste0("3cd2a0c8:", f)), stdout = TRUE),
             file.path(base, f))
}
out <- utils::capture.output(eval(parse(text = src), envir = globalenv()))
cat(out, sep = "\n")
nv <- unique(unlist(strsplit(sub("^[[:space:]]*NEW FREE VARIABLES: ", "",
                                 grep("NEW FREE VARIABLES", out, value = TRUE)),
                             ", ", fixed = TRUE)))
nv <- trimws(nv)
ns <- c(asNamespace("frmtmb"), lapply(c("frmtmb.sample", "frmtmb.ode",
  "frmtmb.survival"), function(p) asNamespace(p)))
for (v in nv) {
  where <- vapply(ns, function(e) exists(v, envir = e, inherits = TRUE), NA)
  cat(sprintf("FREE %-24s found in %s\n", v,
              if (any(where)) paste(c("frmtmb", "sample", "ode", "survival")[where],
                                    collapse = ",") else "NOWHERE"))
}
cat("SUPERASSIGN CHECK DONE\n")
