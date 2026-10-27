# Reviewer of lane loglik, fix round 1: recount the dry-run bin-2 ledger
# against the tracked one. Output: dev/review-loglik-log/fr1-ledger.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-loglik"
rd <- function(f) utils::read.delim(f, quote = "", stringsAsFactors = FALSE,
                                    colClasses = "character")
new <- rd(file.path(wt, "dev/loglik-ledger-out/brmsport-bin2-ledger.tsv"))
old <- rd(file.path(wt, "dev/brmsport-bin2-ledger.tsv"))
b1 <- rd(file.path(wt, "dev/brmsport-ledger.tsv"))
cnt <- function(L, ...) sum(as.integer(L$n[...]))
for (nm in c("old", "new")) {
  L <- get(nm)
  cat(sprintf("%s: rows %d, n total %d, pass 2T %d, pass 2P %d, adapted-input %d\n",
              nm, nrow(L), cnt(L, TRUE),
              cnt(L, L$outcome == "pass" & L$bin == "2T"),
              cnt(L, L$outcome == "pass" & L$bin == "2P"),
              cnt(L, L$outcome == "pass" & L$class == "adapted-input")))
}
cat("bin 1 columns:", paste(head(names(b1), 12), collapse = " "), "\n")
oc <- intersect(c("outcome", "verdict"), names(b1))[1]
nn <- if ("n" %in% names(b1)) as.integer(b1$n) else rep(1L, nrow(b1))
cat(sprintf("bin 1: total %d, pass %d\n", sum(nn), sum(nn[b1[[oc]] == "pass"])))
key <- function(L) paste(L$id, L$file, L$line)
m <- match(key(new), key(old))
cat("new rows not in old:", sum(is.na(m)), "\n")
ch <- which(!is.na(m) & (new$outcome != old$outcome[m] |
                          new$class != old$class[m]))
cat("rows whose outcome or class changed:", length(ch), "\n")
print(data.frame(id = new$id[ch], n = new$n[ch],
                 old = paste(old$outcome[m[ch]], old$class[m[ch]]),
                 new = paste(new$outcome[ch], new$class[ch])),
      row.names = FALSE)
