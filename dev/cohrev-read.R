## The lane's TSV reader, re-implemented for the review so that the
## numbers below are not read back through the code that produced them.
read_kv <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  rows <- lapply(ln, function(l) {
    kv <- strsplit(strsplit(l, "\t", fixed = TRUE)[[1L]], "=",
                   fixed = TRUE)
    setNames(trimws(vapply(kv, function(x) paste(x[-1L], collapse = "="),
                           character(1))),
             trimws(vapply(kv, `[`, character(1), 1L)))
  })
  nms <- unique(unlist(lapply(rows, names)))
  d <- as.data.frame(do.call(rbind, lapply(rows, function(r) r[nms])),
                     stringsAsFactors = FALSE)
  names(d) <- nms
  for (cl in c("rep", "seed", "est", "se", "lo", "hi", "width", "sd_id",
               "sd_idcond", "loglik", "conv", "maxgrad", "nbadse",
               "secs")) {
    if (cl %in% names(d)) {
      d[[cl]] <- suppressWarnings(as.numeric(d[[cl]]))
    }
  }
  for (cl in c("ok", "covers", "pdhess")) {
    if (cl %in% names(d)) d[[cl]] <- d[[cl]] == "TRUE"
  }
  d
}
