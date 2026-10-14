# Lane optstart's exgaussian battery rows (dev/optstart-battery.R, 5
# designs x 10 seeds) on a build of the release tree, against the
# lane's final build (install 13, label lane13) and rellib-r8 (r8): per
# row the log likelihood, the code and the warnings; per design the
# errors and the nonzero codes.
#   Rscript dev/optstart-battery.R <lib> <label> exgaussian
#   Rscript dev/rel071-rr-optstart-exg.R <label>
lab <- commandArgs(TRUE)[1]
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, quote = "")
g <- function(p) {
  fs <- Sys.glob(sprintf("dev/optstart-log/battery/%s-*.tsv", p))
  x <- do.call(rbind, lapply(fs, rd))
  x[x$family == "exgaussian", ]
}
L <- g("lane13"); B <- g("r8")
N <- rd(sprintf("dev/optstart-log/battery/%s-exgaussian.tsv", lab))
M <- merge(L, N, by = c("family", "design", "seed"), suffixes = c(".l", ".n"))
M <- merge(M, B[, c("design", "seed", "ll", "code")], by = c("design", "seed"))
cat("rows: lane", nrow(L), "merged", nrow(N), "matched", nrow(M), "\n")
eqna <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
same_ll <- (is.na(M$ll.l) & is.na(M$ll.n)) |
  (!is.na(M$ll.l) & !is.na(M$ll.n) &
     abs(M$ll.l - M$ll.n) <= 1e-8 * pmax(1, abs(M$ll.l)))
cat("log likelihood the same (1e-8 relative, NA = error on both):",
    sum(same_ll), "of", nrow(M), "\n")
cat("code the same:", sum(eqna(M$code.l, M$code.n)), "of", nrow(M), "\n")
mask <- function(x) gsub("[0-9.e+-]+", "#", ifelse(is.na(x), "", x))
cat("warnings the same (digits masked):",
    sum(mask(M$warn.l) == mask(M$warn.n)), "of", nrow(M), "\n")
tab <- do.call(rbind, lapply(split(seq_len(nrow(M)), M$design), function(i) {
  data.frame(design = M$design[i[1]],
             r8_errors = sum(is.na(M$ll[i])),
             r8_nonzero = sum(M$code[i] != 0, na.rm = TRUE),
             lane_errors = sum(is.na(M$ll.l[i])),
             lane_nonzero = sum(M$code.l[i] != 0, na.rm = TRUE),
             merged_errors = sum(is.na(M$ll.n[i])),
             merged_nonzero = sum(M$code.n[i] != 0, na.rm = TRUE))
}))
print(tab, row.names = FALSE)
cat("error seeds, rebig: lane", M$seed[M$design == "rebig" & is.na(M$ll.l)],
    "| merged", M$seed[M$design == "rebig" & is.na(M$ll.n)], "\n")
d <- M[which(!same_ll | !eqna(M$code.l, M$code.n)),
       c("design", "seed", "ll.l", "ll.n", "code.l", "code.n")]
cat("rows that differ from the lane:", nrow(d), "\n")
if (nrow(d)) print(d, row.names = FALSE)
