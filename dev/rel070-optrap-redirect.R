# Release 0.70.0, rerun of lane optrap's studies: run one lane or
# reviewer script unchanged, with the files it writes moved into
# dev/rel070-log/studies/optrap/, because the release tree must not get
# new files in the lane's and reviewer's log directories. Modeled on
# dev/optrap-p2-redirect.R.
#   Rscript dev/rel070-optrap-redirect.R <script> [its args ...]
args0 <- commandArgs(TRUE)
out_root <- "dev/rel070-log/studies/optrap"
map_path <- function(p) {
  if (!is.character(p)) return(p)
  p <- sub("^dev/optrap-rev3-log", file.path(out_root, "rev3-log"), p)
  p <- sub("^dev/optrap-rev2-log", file.path(out_root, "rev2-log"), p)
  p <- sub("^dev/optrap-rev-log", file.path(out_root, "rev-log"), p)
  sub("^dev/optrap-log", file.path(out_root, "optrap-log"), p)
}
for (d in c("rev3-log", "rev2-log", "rev-log", "optrap-log")) {
  base::dir.create(file.path(out_root, d), recursive = TRUE,
                   showWarnings = FALSE)
}
write.table <- function(x, file = "", ...) {
  utils::write.table(x, map_path(file), ...)
}
dir.create <- function(path, ...) base::dir.create(map_path(path), ...)
a <- vapply(args0[-1], map_path, "", USE.NAMES = FALSE)
commandArgs <- function(trailingOnly = FALSE) a
source(args0[1])
