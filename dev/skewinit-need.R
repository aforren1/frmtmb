# Which packages the lane still needs, and which of them are hollow.
U <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
LANE <- "C:/Users/adf44/source/r/skewinit-lib"
REL <- "C:/Users/adf44/source/r/rellib-r3"
have <- function(L) {
  if (!dir.exists(L)) return(character(0))
  d <- list.dirs(L, recursive = FALSE, full.names = FALSE)
  d[file.exists(file.path(L, d, "DESCRIPTION"))]
}
ok <- unique(c(have(LANE), have(REL), have(U),
               rownames(utils::installed.packages(
                 lib.loc = .Library))))
dcf <- read.dcf(file.path(REL, "frmtmb", "DESCRIPTION"))
split_dep <- function(f) {
  if (!f %in% colnames(dcf)) return(character(0))
  x <- strsplit(dcf[1, f], ",")[[1]]
  x <- trimws(sub("\\(.*", "", x))
  setdiff(x, c("R", ""))
}
imp <- unique(c(split_dep("Imports"), split_dep("Depends"),
                split_dep("LinkingTo")))
sug <- split_dep("Suggests")
cat("frmtmb Imports/Depends/LinkingTo:", length(imp), "\n")
cat("  missing:", paste(setdiff(imp, ok), collapse = " "), "\n")
cat("frmtmb Suggests:", length(sug), "\n")
cat("  missing:", paste(setdiff(sug, ok), collapse = " "), "\n")
cat("\nhollow in the user library:",
    length(setdiff(list.dirs(U, recursive = FALSE, full.names = FALSE),
                   have(U))), "\n")
