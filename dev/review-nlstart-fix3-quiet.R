# Review of lane nlstart, third short check: fresh rows where a probe win
# returns a short fit with no warning while 0.79.1 (= the first fit) warned.
a <- commandArgs(trailingOnly = TRUE)
rd1 <- function(dir) {
  fs <- list.files(dir, "[.]log$", full.names = TRUE)
  ln <- unlist(lapply(fs, function(f) grep("^ROW ", readLines(f), value = TRUE)))
  p <- strsplit(ln, " ")
  data.frame(fam = sapply(p, `[`, 2), body = sapply(p, `[`, 3),
             seed = sapply(p, `[`, 4),
             ll = suppressWarnings(as.numeric(sapply(p, `[`, 5))),
             nwarn = as.integer(sapply(p, `[`, 8)), wclass = sapply(p, `[`, 9),
             multi = as.integer(sapply(p, `[`, 10)),
             ref = suppressWarnings(as.numeric(sapply(p, `[`, 14))))
}
rd <- function(dirs) do.call(rbind, lapply(seq_along(dirs), function(k) {
  d <- rd1(dirs[k]); d$set <- k; d }))
A <- rd(strsplit(a[1], ",")[[1]]); B <- rd(strsplit(a[2], ",")[[1]])
key <- function(d) paste(d$set, d$fam, d$body, d$seed)
B <- B[match(key(A), key(B)), ]
best <- pmax(A$ref, A$ll, B$ll, na.rm = TRUE)
sh <- !is.na(A$ll) & best - A$ll > 1e-4
q <- sh & A$nwarn == 0 & A$multi == 1 & B$nwarn > 0
cat(sprintf("rows %d; probe wins %d; short after a probe win %d; of them silent %d, silent where 0.79.1 warned %d\n",
            nrow(A), sum(A$multi == 1, na.rm = TRUE), sum(sh & A$multi == 1),
            sum(sh & A$multi == 1 & A$nwarn == 0), sum(q)))
print(data.frame(A[q, c("fam", "body", "seed", "ll")], base = B$ll[q],
                 base_nwarn = B$nwarn[q], short = (best - A$ll)[q]), row.names = FALSE)
cat("warning classes on the lane rows:", paste(sort(unique(A$wclass)), collapse = " "), "\n")
