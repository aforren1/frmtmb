# Restore the R user library after a file-level sweep. See
# dev/machine-library.md for the diagnosis and the standing decisions.
LIB <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(LIB)
options(repos = c(CRAN = "https://cloud.r-project.org"))
dirs <- list.dirs(LIB, recursive = FALSE, full.names = FALSE)
hollow <- dirs[!file.exists(file.path(LIB, dirs, "DESCRIPTION"))]
hollow <- setdiff(hollow, grep("^(00LOCK|frmtmb)", hollow, value = TRUE))
av <- rownames(available.packages(type = "binary"))
todo <- intersect(hollow, av)
cat("hollow:", length(hollow), " on CRAN binary:", length(todo),
    " not on CRAN:", length(setdiff(hollow, av)), "\n")
cat("NOT ON CRAN:", paste(setdiff(hollow, av), collapse = " "), "\n")
if (length(todo)) {
  install.packages(todo, lib = LIB, type = "binary", dependencies = FALSE)
}
after <- file.exists(file.path(LIB, todo, "DESCRIPTION"))
cat("RESTORED", sum(after), "of", length(todo), "\n")
cat("STILL HOLLOW:", paste(todo[!after], collapse = " "), "\n")
cat("RESTORE DONE\n")
