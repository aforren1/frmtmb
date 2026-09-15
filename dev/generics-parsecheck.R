# After an interrupted session: does every R file this lane touched
# still parse, is every test file whole, and is the lane's library
# free of hollow package directories?  Run before anything else,
# because a process killed mid-edit can leave a source half-written.
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
changed <- system2("git", c("-C", shQuote(root), "status", "--porcelain"),
                   stdout = TRUE)
paths <- sub("^...", "", changed)
paths <- paths[grepl("[.]R$", paths)]
extra <- c(list.files(file.path(root, "R"), full.names = FALSE,
                      pattern = "[.]R$"))
paths <- unique(c(paths, file.path("R", extra),
                  "tests/testthat/test-generic-collision.R",
                  "tests/testthat/test-scale-contract.R"))
paths <- paths[!grepl("/$", paths)]
bad <- character()
n <- 0L
for (p in paths) {
  f <- file.path(root, p)
  if (dir.exists(f) || !file.exists(f)) next
  n <- n + 1L
  r <- tryCatch({ parse(f); "ok" }, error = function(e) conditionMessage(e))
  if (!identical(r, "ok")) bad <- c(bad, paste(p, r))
}
cat(sprintf("R files parsed: %d, failing: %d\n", n, length(bad)))
for (b in bad) cat("  ", b, "\n")
# the untracked dev scripts too
dv <- list.files(file.path(root, "dev"), pattern = "^generics-.*[.]R$",
                 full.names = TRUE)
db <- Filter(function(f) !identical(tryCatch({ parse(f); "ok" },
  error = function(e) "x"), "ok"), dv)
cat(sprintf("dev/generics-*.R parsed: %d, failing: %d\n", length(dv),
            length(db)))
for (b in db) cat("  ", basename(b), "\n")
# a hollow package directory is one with no DESCRIPTION or no R/ db
for (lib in c("C:/Users/adf44/source/r/generics-lib")) {
  d <- list.dirs(lib, recursive = FALSE)
  hollow <- d[!file.exists(file.path(d, "DESCRIPTION")) |
                !file.exists(file.path(d, "R", basename(d)))]
  cat(sprintf("%s: %d packages, %d hollow\n", lib, length(d),
              length(hollow)))
  for (h in hollow) cat("  ", h, "\n")
}
