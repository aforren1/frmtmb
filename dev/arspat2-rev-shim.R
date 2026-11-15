# Fix round: the reviewer's scripts read autocor_seam_t_missing(), which
# the fix removed. This stand-in reports the margin route the way they
# expect: a registered term's response whose missing values left the fit.
seam_tm <- function(frame) {
  out <- list()
  for (r in names(frame[["mi_map"]])) {
    ac <- frame[["autocor"]][[r]]
    un <- frame[["mi_map"]][[r]][["unread"]]
    if (!is.null(ac) && isTRUE(ac[["seam"]]) && length(un)) {
      out[[r]] <- list(m = length(un), rows = un, margin = TRUE)
    }
  }
  out
}
