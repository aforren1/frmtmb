# Reviewer scratch: diff the two dumps of dev/review-compat-dump.R.
b <- readRDS("dev/review-compat-before.rds")
a <- readRDS("dev/review-compat-after.rds")
key <- function(x) paste(pmin(x$feature_a, x$feature_b),
                         pmax(x$feature_a, x$feature_b), sep = " x ")
for (w in c("core", "all")) {
  B <- b[[w]]; A <- a[[w]]
  B$key <- key(B); A$key <- key(A)
  m <- merge(B[c("key", "status", "note")], A[c("key", "status", "note")],
             by = "key", suffixes = c(".b", ".a"))
  cat("\n==", w, ": pairs", nrow(B), nrow(A), "merged", nrow(m), "\n")
  st <- m[m$status.b != m$status.a, ]
  nt <- m[m$status.b == m$status.a & m$note.b != m$note.a, ]
  cat("status moves", nrow(st), "; note-only", nrow(nt), "\n")
  print(table(paste(st$status.b, "->", st$status.a)))
  if (w == "all") {
    cat("moves outside core keys:\n")
    ck <- key(a$core)
    print(st[!st$key %in% ck, c("key", "status.b", "status.a")])
    print(nt[!nt$key %in% ck, c("key", "status.b", "status.a")])
  }
}
A <- a$all
ft <- a$ft
ats <- ft$name[ft$kind == "aterm"]
only <- c("thres()", "disc", "equidistant", "sum_to_zero", "cs_pred()")
cat("\naterm x ordinal_only, all loaded, after:\n")
for (at in ats) {
  s <- vapply(only, function(o) {
    r <- A$status[(A$feature_a == at & A$feature_b == o) |
                    (A$feature_a == o & A$feature_b == at)]
    if (length(r)) r[1] else "-"
  }, "")
  cat(sprintf("%-10s %s\n", at, paste(sprintf("%s=%s", only, s),
                                      collapse = " ")))
}
cat("\naterm kinds and sources:\n")
print(ft[ft$kind == "aterm", ])
