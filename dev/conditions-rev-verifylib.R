# Reviewer, lane wt-conditions: does the lane build match the worktree
# source? For every top-level `name <- function` in each package's R/,
# compare the parsed function with the installed namespace object.
LIB <- Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/conditions-lib")
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
pkgs <- c(frmtmb = root,
          setNames(file.path(root, "extensions", list.files(
            file.path(root, "extensions"))),
            list.files(file.path(root, "extensions"))))
norm <- function(f) {
  f <- removeSource(f)
  environment(f) <- globalenv()
  paste(deparse(f, control = c("keepInteger", "niceNames")),
        collapse = "\n")
}
for (p in names(pkgs)) {
  ns <- asNamespace(p)
  cat(p, "installed from", find.package(p), "\n")
  files <- list.files(file.path(pkgs[[p]], "R"), full.names = TRUE)
  n <- 0; bad <- character(); missing <- character()
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.name(e[[2]]) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) {
        nm <- as.character(e[[2]])
        src <- eval(e[[3]], globalenv())
        if (!exists(nm, envir = ns, inherits = FALSE)) {
          missing <- c(missing, nm); next
        }
        obj <- get(nm, envir = ns, inherits = FALSE)
        if (!is.function(obj)) next
        n <- n + 1
        if (!identical(norm(src), norm(obj))) bad <- c(bad, paste(basename(f), nm))
      }
    }
  }
  cat(sprintf("  functions compared %d, differing %d, not in ns %d\n",
              n, length(bad), length(missing)))
  if (length(bad)) print(head(bad, 20))
}
