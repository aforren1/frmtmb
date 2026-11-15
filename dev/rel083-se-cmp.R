# Lane selimit2's SE comparison of the merged build (rellib-r22) against
# 0.82.0 (rellib-r21), by record identity: a record is named by its file,
# its test's description, its kind (COV, VCT) and its rank among the
# records of that kind in that test. Records only one arm has are counted,
# not compared. Then every firing of the limit rule on the merged build
# (a LIM line, which belongs to the COV record after it), classified as
# lane selimit2's dev/boundfam-classify.R does: the SEs it NaNs were NaN
# on 0.82.0 already (already_nan), finite there (newly_nan), or the
# record has no 0.82.0 twin (no_base).
#
#   Rscript dev/rel083-se-cmp.R dev/rel083-log/se-r22 dev/rel083-log/se-r21
a <- commandArgs(trailingOnly = TRUE)
dn <- a[1]; db <- a[2]
num <- function(s) suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1]]))
keyed <- function(f) {
  if (!file.exists(f)) return(list())
  x <- readLines(f, warn = FALSE)
  test <- "(file level)"; seen <- list(); out <- list(); lim <- NULL
  for (l in x) {
    if (startsWith(l, "TEST ")) {
      test <- sub("^TEST ", "", l); seen <- list(); next
    }
    tag <- sub(" .*", "", l)
    if (tag == "LIM") { lim <- l; next }
    if (!tag %in% c("COV", "VCT")) next
    seen[[tag]] <- (seen[[tag]] %||% 0L) + 1L
    k <- paste(test, tag, seen[[tag]], sep = " | ")
    v <- if (tag == "COV") sub(".*se=", "", l) else sub(".*se_t=", "", l)
    out[[k]] <- list(v = num(v), tag = tag, lim = if (tag == "COV") lim,
                     lost = if (tag == "COV") sub(" se=.*", "", sub(".*lost=", "", l)))
    if (tag == "COV") lim <- NULL
  }
  out
}
`%||%` <- function(x, y) if (is.null(x)) y else x
cmp1 <- function(x, y) {
  if (length(x) != length(y)) return("length")
  if (!identical(is.na(x), is.na(y))) return("nan")
  if (identical(x, y)) return("bitwise")
  ok <- !is.na(x)
  r <- max(abs(x[ok] - y[ok]) / pmax(abs(x[ok]), 1e-300))
  if (r <= 1e-8) "ulps" else "moved"
}
files <- basename(Sys.glob(file.path(dn, "*.trace")))
libs <- vapply(c(dn, db), function(d) {
  l <- unlist(lapply(Sys.glob(file.path(d, "*.trace")), readLines, n = 1L))
  paste(names(table(l)), table(l), collapse = "; ")
}, "")
cat("trace lib lines:\n "); cat(libs, sep = "\n ")
lv <- c("bitwise", "ulps", "moved", "nan", "length", "only_new", "only_base")
tot <- list(COV = character(0), VCT = character(0))
cls <- c(fired = 0L, already_nan = 0L, newly_nan = 0L, no_base = 0L)
for (f in files) {
  n <- keyed(file.path(dn, f)); b <- keyed(file.path(db, f))
  ks <- union(names(n), names(b))
  res <- vapply(ks, function(k) {
    if (is.null(b[[k]])) return("only_new")
    if (is.null(n[[k]])) return("only_base")
    cmp1(n[[k]]$v, b[[k]]$v)
  }, "")
  tag <- vapply(ks, function(k) (n[[k]] %||% b[[k]])$tag, "")
  for (t in names(tot)) tot[[t]] <- c(tot[[t]], res[tag == t])
  bad <- ks[res %in% c("moved", "nan", "length")]
  if (length(bad)) {
    cat(sprintf("%s: %d of %d keyed records differ\n", f, length(bad), length(ks)))
    for (k in head(bad, 4)) {
      cat("   ", k, "[", res[[k]], "]\n     new: ", signif(n[[k]]$v, 4),
          "\n     0.82.0:", signif(b[[k]]$v, 4), "\n")
    }
  }
  for (k in names(n)) {
    if (is.null(n[[k]]$lim)) next
    cls["fired"] <- cls["fired"] + 1L
    cat("LIM ", f, " | ", k, "\n   ", substr(n[[k]]$lim, 1, 300), "\n", sep = "")
    bb <- b[[k]]
    if (is.null(bb) || length(bb$v) != length(n[[k]]$v)) {
      cls["no_base"] <- cls["no_base"] + 1L; cat("   no base\n"); next
    }
    new <- which(is.na(n[[k]]$v) & !is.na(bb$v))
    if (length(new)) {
      cls["newly_nan"] <- cls["newly_nan"] + 1L
      cat("   0.82.0 SE -> NaN at", new, ":", signif(bb$v[new], 4), "| lost", n[[k]]$lost, "\n")
    } else cls["already_nan"] <- cls["already_nan"] + 1L
  }
}
cat("files", length(files), "\n")
for (t in names(tot)) {
  tb <- table(factor(tot[[t]], levels = lv))
  cat(t, paste(names(tb), tb, sep = "=", collapse = " "), "\n")
}
cat("TOTAL", paste(names(cls), cls, sep = "=", collapse = " "), "\n")
