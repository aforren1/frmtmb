# Review of lane readers: compare regress-before.rds and regress-after.rds.
# "same" is identical(); otherwise the largest relative difference, and for
# conditional_effects() each column separately.
a <- readRDS("dev/review-readers-out/regress-before.rds")
b <- readRDS("dev/review-readers-out/regress-after.rds")
rdiff <- function(x, y) {
  if (identical(x, y)) return("same")
  if (is.character(x) || is.character(y)) {
    return(paste0("before: ", substr(paste(x, collapse = ""), 1, 70),
                  " | after: ", substr(paste(y, collapse = ""), 1, 70)))
  }
  ux <- unlist(x); uy <- unlist(y)
  if (length(ux) != length(uy)) return("length differs")
  num <- is.numeric(ux) && is.numeric(uy)
  if (!num) return("non-numeric differs")
  nax <- is.na(ux); nay <- is.na(uy)
  if (any(nax != nay)) {
    return(sprintf("NA pattern differs (before NA %d, after NA %d)",
                   sum(nax), sum(nay)))
  }
  k <- !nax
  r <- max(abs(uy[k] - ux[k]) / pmax(abs(ux[k]), 1e-300))
  sprintf("max rel %.2e", r)
}
for (nm in names(b)) {
  for (w in setdiff(names(b[[nm]]), "ce")) {
    cat(sprintf("%-16s %-11s %s\n", nm, w, rdiff(a[[nm]][[w]], b[[nm]][[w]])))
  }
  ca <- a[[nm]]$ce; cb <- b[[nm]]$ce
  if (is.character(ca) || is.character(cb)) {
    cat(sprintf("%-16s %-11s %s\n", nm, "ce", rdiff(ca, cb)))
    next
  }
  for (k in names(cb)) {
    for (col in intersect(c("estimate__", "se__", "lower__", "upper__"),
                          names(cb[[k]]))) {
      cat(sprintf("%-16s ce %-8s %-10s %s\n", nm, k, col,
                  rdiff(ca[[k]][[col]], cb[[k]][[col]])))
    }
  }
}
