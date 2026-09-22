# Reviewer, priority 3: every CALL to predict(), fixef(), vcov(),
# fitted(), residuals() and ngrps() in every R file of all eight
# packages, found by PARSING rather than by grep, so that a mention
# inside a string or a comment cannot be counted and a call split over
# two lines cannot be missed.
#
#   Rscript dev/shapes-rev-callsites.R > dev/shapes-rev-callsites.txt

TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
TARGETS <- c("predict", "fixef", "vcov", "fitted", "residuals", "ngrps",
             "nsamples", "posterior_samples", "parnames")

walk <- function(e, file, out) {
  if (is.call(e)) {
    fn <- e[[1L]]
    nm <- NULL
    if (is.name(fn)) nm <- as.character(fn)
    else if (is.call(fn) && length(fn) == 3L &&
             as.character(fn[[1L]]) %in% c("::", ":::"))
      nm <- as.character(fn[[3L]])
    if (!is.null(nm) && nm %in% TARGETS) {
      txt <- paste(deparse(e), collapse = " ")
      out$hits[[length(out$hits) + 1L]] <-
        list(file = file, fn = nm, txt = substr(txt, 1, 160))
    }
    for (i in seq_along(e)) {
      if (!is.null(e[[i]]) && !identical(e[[i]], quote(expr = )))
        out <- walk(e[[i]], file, out)
    }
  } else if (is.pairlist(e) || is.expression(e) ||
             (is.list(e) && !is.function(e))) {
    for (i in seq_along(e)) out <- walk(e[[i]], file, out)
  }
  out
}

dirs <- c(file.path(TREE, "R"),
          Sys.glob(file.path(TREE, "extensions/*/R")))
files <- unlist(lapply(dirs, list.files, pattern = "[.][Rr]$",
                       full.names = TRUE))
out <- list(hits = list())
for (f in files) {
  ex <- tryCatch(parse(f, keep.source = FALSE),
                 error = function(e) { cat("PARSE FAIL ", f, "\n"); NULL })
  if (is.null(ex)) next
  for (i in seq_along(ex)) out <- walk(ex[[i]], f, out)
}
df <- do.call(rbind, lapply(out$hits, function(h)
  data.frame(file = sub(paste0(TREE, "/"), "", h$file, fixed = TRUE),
             fn = h$fn, txt = h$txt, stringsAsFactors = FALSE)))
cat("total call sites:", nrow(df), "\n\n")
for (fn in TARGETS) {
  s <- df[df$fn == fn, ]
  if (!nrow(s)) next
  cat("======== ", fn, " (", nrow(s), ") ========\n", sep = "")
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-52s %s\n", s$file[i], s$txt[i]))
  cat("\n")
}
