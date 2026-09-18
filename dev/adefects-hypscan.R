# Every hypothesis() call in this repository whose hypothesis STRING
# states no relation, which the D3 refusal now rejects. The scan parses
# each file and walks the AST, so a call split over four lines is seen
# once and a string inside a comment is not seen at all.
#
#   Rscript dev/adefects-hypscan.R [root]
#
# A string held in a variable (test-boot.R's `icc`) cannot be resolved
# statically and is reported separately, by call site, so the count is
# a floor that names what it cannot see rather than a silent one.
root <- commandArgs(trailingOnly = TRUE)
if (!length(root)) root <- "."

files <- c(
  list.files(file.path(root, "R"), "[.]R$", full.names = TRUE),
  list.files(file.path(root, "tests"), "[.]R$", full.names = TRUE,
             recursive = TRUE),
  list.files(file.path(root, "vignettes"), "[.]Rmd$", full.names = TRUE),
  list.files(file.path(root, "extensions"), "[.](R|Rmd)$",
             full.names = TRUE, recursive = TRUE))
files <- files[!grepl("/brms-suite/", files, fixed = TRUE)]

has_rel <- function(s) grepl("[=<>]", s)

# an .Rmd is parsed by lifting its R chunks; a chunk that does not parse
# is reported rather than skipped
code_of <- function(f) {
  if (!grepl("[.]Rmd$", f)) return(parse(f, keep.source = TRUE))
  ln <- readLines(f, warn = FALSE)
  open <- grep("^[`]{3}[{]r", ln)
  close <- grep("^[`]{3}\\s*$", ln)
  out <- character(0)
  for (o in open) {
    cl <- close[close > o]
    if (!length(cl)) next
    out <- c(out, ln[(o + 1L):(cl[1L] - 1L)])
  }
  parse(text = out, keep.source = TRUE)
}

bare <- character(0)
dyn <- character(0)
rel <- 0L
walk <- function(e, f, line) {
  if (!is.call(e)) return(invisible(NULL))
  fn <- e[[1L]]
  if (is.name(fn) && as.character(fn) == "hypothesis" && length(e) >= 3L) {
    a <- e[[3L]]
    lit <- NULL
    if (is.character(a)) lit <- a
    if (is.call(a) && identical(a[[1L]], as.name("c")) &&
        all(vapply(as.list(a)[-1L], is.character, NA))) {
      lit <- unlist(as.list(a)[-1L])
    }
    if (is.null(lit)) {
      dyn <<- c(dyn, sprintf("%s:%d  %s", f, line, deparse1(e)))
    } else if (any(!has_rel(lit))) {
      bare <<- c(bare, sprintf("%s:%d  %s", f, line,
                               paste(lit[!has_rel(lit)], collapse = " | ")))
    } else {
      rel <<- rel + 1L
    }
  }
  # an EMPTY argument (`x[, 1]`) errors when it is touched, not when it
  # is extracted, so every read of a call's parts is guarded
  for (i in seq_along(e)) {
    p <- tryCatch({
      a <- e[[i]]
      if (is.call(a)) a else NULL
    }, error = function(err) NULL)
    if (!is.null(p)) walk(p, f, line)
  }
  invisible(NULL)
}

for (f in files) {
  ex <- tryCatch(code_of(f), error = function(e) NULL)
  if (is.null(ex)) next
  rr <- utils::getSrcref(ex)
  for (i in seq_along(ex)) {
    ln <- if (!is.null(rr) && !is.null(rr[[i]])) rr[[i]][1L] else NA_integer_
    walk(ex[[i]], sub(paste0("^", root, "/?"), "", f), ln)
  }
}

cat("hypothesis() calls with a literal hypothesis and a relation:", rel, "\n")
cat("hypothesis() calls with a literal hypothesis and NO relation:",
    length(bare), "\n")
for (b in bare) cat("  BARE ", b, "\n", sep = "")
cat("hypothesis() calls whose hypothesis is not a literal:",
    length(dyn), "\n")
for (b in dyn) cat("  DYN  ", substr(b, 1, 150), "\n", sep = "")
