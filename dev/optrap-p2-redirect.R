# Lane optrap, punch round 2: run one of the reviewer's scripts with its
# outputs redirected into this lane's dev/optrap-log/ (the reviewer's
# dev/optrap-rev-log/ and dev/optrap-rev2-log/ are theirs). Paths under
# those directories given to write.table(), read.delim() (for the
# reviewer's own intermediate files), dir.create() and file arguments
# are mapped to dev/optrap-log/rev/ and dev/optrap-log/rev2/.
#   Rscript dev/optrap-p2-redirect.R <script> [its args ...]
args0 <- commandArgs(TRUE)
map_path <- function(p) {
  if (!is.character(p)) return(p)
  p <- sub("^dev/optrap-rev2-log", "dev/optrap-log/rev2", p)
  sub("^dev/optrap-rev-log", "dev/optrap-log/rev", p)
}
dir.create("dev/optrap-log/rev2", recursive = TRUE, showWarnings = FALSE)
dir.create("dev/optrap-log/rev", recursive = TRUE, showWarnings = FALSE)
write.table <- function(x, file = "", ...) {
  utils::write.table(x, map_path(file), ...)
}
dir.create <- function(path, ...) base::dir.create(map_path(path), ...)
a <- vapply(args0[-1], map_path, "", USE.NAMES = FALSE)
commandArgs <- function(trailingOnly = FALSE) a
source(args0[1])
