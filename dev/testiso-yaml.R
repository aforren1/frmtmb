# Parse .github/workflows/test-coverage.yaml, list its jobs and matrix,
# check every Rscript step parses as R, and check each matrix entry's
# siblings against the frmtmb packages its DESCRIPTION names.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
y <- yaml::read_yaml(".github/workflows/test-coverage.yaml")
cat("jobs:", names(y$jobs), "\n")
for (j in names(y$jobs)) {
  for (s in y$jobs[[j]]$steps) {
    if (identical(s$shell, "Rscript {0}")) {
      ok <- tryCatch({
        parse(text = s$run)
        "parses"
      }, error = function(e) conditionMessage(e))
      cat(j, "/", s$name, ":", ok, "\n")
    }
  }
  cat(j, "cache-version:", unlist(lapply(y$jobs[[j]]$steps, function(s) {
    s$with$`cache-version`
  })), "\n")
}
named <- function(pkg) {
  d <- read.dcf(file.path("extensions", pkg, "DESCRIPTION"))
  f <- intersect(c("Depends", "Imports", "Suggests"), colnames(d))
  x <- trimws(sub("[(].*", "", unlist(strsplit(d[, f], ","))))
  setdiff(x[startsWith(x, "frmtmb")], "frmtmb")
}
imports <- function(pkg) {
  d <- read.dcf(file.path("extensions", pkg, "DESCRIPTION"))
  f <- intersect(c("Depends", "Imports"), colnames(d))
  x <- trimws(sub("[(].*", "", unlist(strsplit(d[, f], ","))))
  setdiff(x[startsWith(x, "frmtmb")], "frmtmb")
}
inc <- y$jobs[["extension-coverage"]]$strategy$matrix$include
cat("matrix:", vapply(inc, `[[`, "", "pkg"), "\n")
stopifnot(setequal(vapply(inc, `[[`, "", "pkg"), list.files("extensions")))
for (e in inc) {
  sib <- strsplit(trimws(e$siblings), " +")[[1]]
  sib <- sib[nzchar(sib)]
  # every sibling's own hard frmtmb dependencies come before it
  order_ok <- all(vapply(seq_along(sib), function(i) {
    all(imports(sib[i]) %in% sib[seq_len(i - 1L)])
  }, TRUE))
  cat(e$pkg, "| names:", named(e$pkg), "| siblings:", sib,
      "| same set:", setequal(named(e$pkg), sib),
      "| order ok:", order_ok, "\n")
}
# the variables the matrix job's tests read, set as the check jobs set them
ev <- y$jobs[["extension-coverage"]]$env
cat("extension-coverage env:", paste0(names(ev), "=", unlist(ev)), "\n")
stopifnot(identical(ev$FRMTMB_SAMPLER_GATES, "false"),
          identical(ev$NOT_CRAN, "true"))
