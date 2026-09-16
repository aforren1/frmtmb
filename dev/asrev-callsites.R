## REVIEW: frm_check_dots() reads the CALLING FRAME. Is every call site
## the method itself, or does any of them sit in a helper the method
## calls? A helper frame reports the helper's name and the helper's
## formals, which is a wrong message.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-argspell"
dirs <- c(file.path(wt, "R"), file.path(wt, "extensions/frmtmb.sample/R"))
tot <- 0L
suspect <- character()
for (d in dirs) {
  for (fp in list.files(d, pattern = "[.]R$", full.names = TRUE)) {
    ex <- parse(fp)
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      if (!is.call(e) ||
          !as.character(e[[1L]])[1L] %in% c("<-", "=")) next
      rhs <- e[[3L]]
      if (!is.call(rhs) ||
          !identical(as.character(rhs[[1L]]), "function")) next
      nm <- as.character(e[[2L]])
      b <- rhs[[3L]]
      txt <- deparse(b)
      hits <- grep("frm_check_dots", txt, fixed = TRUE)
      if (!length(hits)) next
      tot <- tot + length(hits)
      # the call must be a TOP-LEVEL statement of the function body,
      # not nested inside another function defined within it
      top <- FALSE
      if (is.call(b) && identical(as.character(b[[1L]]), "{")) {
        for (j in seq_along(b)[-1L]) {
          s <- b[[j]]
          if (is.call(s) &&
              identical(as.character(s[[1L]])[1L], "frm_check_dots")) {
            top <- TRUE
          }
        }
      } else if (is.call(b) &&
                 identical(as.character(b[[1L]])[1L], "frm_check_dots")) {
        top <- TRUE
      }
      has_dots <- "..." %in% names(rhs[[2L]])
      if (!top || !has_dots) {
        suspect <- c(suspect, sprintf("%s  %s  top=%s dots=%s",
                                      basename(fp), nm, top, has_dots))
      }
    }
  }
}
cat("frm_check_dots call sites found in a top-level assignment:", tot, "\n")
cat("sites NOT a top-level statement of a dots-taking function:",
    length(suspect), "\n")
if (length(suspect)) cat(paste0("  ", suspect), sep = "\n")
