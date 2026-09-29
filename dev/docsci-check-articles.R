# Every name in an `articles:` block of every _pkgdown.yml, against the
# .Rmd files that exist. A name with no file would be an article pkgdown
# cannot build; a file with no name would be an article the navbar does
# not list.
dirs <- c(".", Sys.glob("extensions/*"))
dirs <- dirs[file.exists(file.path(dirs, "_pkgdown.yml"))]
for (d in dirs) {
  cfg <- yaml::read_yaml(file.path(d, "_pkgdown.yml"))
  named <- unlist(lapply(cfg$articles, function(s) s$contents))
  named <- named[nzchar(named)]
  files <- sub("[.]Rmd$", "", basename(Sys.glob(
    file.path(d, "vignettes", "*.Rmd"))))
  cat("== ", d, "\n", sep = "")
  cat("   named in yml: ", length(named), "  .Rmd on disk: ",
      length(files), "\n", sep = "")
  m1 <- setdiff(named, files)
  m2 <- setdiff(files, named)
  if (length(m1)) cat("   NAMED WITH NO FILE: ",
                      paste(m1, collapse = " "), "\n", sep = "")
  if (length(m2)) cat("   FILE NOT NAMED: ",
                      paste(m2, collapse = " "), "\n", sep = "")
  if (!length(m1) && !length(m2)) cat("   match\n")
}
