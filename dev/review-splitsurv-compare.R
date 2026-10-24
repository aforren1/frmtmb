# Compare the before (rellib-r11) and after (lane build) dumps of
# dev/splitsurv-curve-dump.R and dev/splitsurv-compat-dump.R.
#   Rscript dev/splitsurv-compare.R
O <- "dev/review-splitsurv-out"
cat("== curve outputs\n")
b <- readRDS(file.path(O, "curve-before.rds"))
a <- readRDS(file.path(O, "curve-after.rds"))
stopifnot(identical(names(a), names(b)))
ulps <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(0)
  max(abs(x[ok] - y[ok]) / pmax(abs(y[ok]), .Machine$double.xmin)) /
    .Machine$double.eps
}
for (nm in names(b)) {
  x <- a[[nm]]; y <- b[[nm]]
  if (inherits(y, "dump_error")) {
    cat(sprintf("%-18s error message identical: %s\n", nm,
                identical(unclass(x), unclass(y))))
    next
  }
  if (is.data.frame(y)) {
    cols <- identical(as.data.frame(unclass(x))[names(y)],
                      as.data.frame(unclass(y)))
    attrs <- setdiff(names(attributes(y)), c("fit", "names", "row.names",
                                             "class"))
    same_attr <- vapply(attrs, function(k) {
      identical(attr(x, k), attr(y, k))
    }, NA)
    rel <- max(vapply(names(y), function(k) {
      if (is.numeric(y[[k]])) ulps(x[[k]], y[[k]]) else 0
    }, 0))
    cat(sprintf("%-18s columns identical: %s; attributes identical: %s; %s\n",
                nm, cols, all(same_attr),
                if (!cols || !all(same_attr)) {
                  paste0("max rel diff in ulps ", format(rel, digits = 3),
                         "; differing attrs: ",
                         paste(attrs[!same_attr], collapse = ","))
                } else "bitwise"))
  } else {
    cat(sprintf("%-18s identical: %s\n", nm, identical(x, y)))
  }
}
cat("== frm_compat()\n")
cb <- readRDS(file.path(O, "compat-before.rds"))
ca <- readRDS(file.path(O, "compat-after.rds"))
key_tab <- function(t) {
  t <- t[order(t$key), c("key", "status", "note")]
  rownames(t) <- NULL
  t
}
tb <- key_tab(cb$tab); ta <- key_tab(ca$tab)
cat("pairs before", nrow(tb), "after", nrow(ta), "; same pair set:",
    identical(tb$key, ta$key), "\n")
watch <- c("cox", "royston_parmar", "rp_floored", "frm_curve",
           "frm_curve_contrast")
m <- merge(tb, ta, by = "key", all = TRUE, suffixes = c(".b", ".a"))
m$sides <- strsplit(m$key, " x ", fixed = TRUE)
inw <- vapply(m$sides, function(s) any(s %in% watch), NA)
diff_status <- m$status.b != m$status.a | is.na(m$status.b) |
  is.na(m$status.a)
diff_note <- !diff_status & m$note.b != m$note.a
for (w in watch) {
  sel <- vapply(m$sides, function(s) w %in% s, NA)
  cat(sprintf("%-20s pairs %d; status differs %d; note differs %d\n", w,
              sum(sel), sum(diff_status & sel), sum(diff_note & sel)))
}
cat("all pairs: status differs", sum(diff_status), "; note differs",
    sum(diff_note), "\n")
if (any(diff_status | diff_note)) {
  d <- m[diff_status | diff_note, c("key", "status.b", "status.a")]
  print(d, row.names = FALSE)
}
fb <- cb$features[order(cb$features$name), c("name", "kind")]
fa <- ca$features[order(ca$features$name), c("name", "kind")]
rownames(fb) <- rownames(fa) <- NULL
cat("features identical (name, kind):", identical(fb, fa), "\n")
cat("unresolved before:", length(cb$unresolved), "after:",
    length(ca$unresolved), "\n")
