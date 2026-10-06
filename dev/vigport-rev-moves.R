# Reviewer: every expression whose class got worse against 0.34.0,
# including "other" rows and rows that ran only patched then.
#
#   Rscript dev/vigport-rev-moves.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
M <- readRDS(file.path(root, "vigport-port-out/r5/results-merged.rds"))
B <- readRDS(file.path(root, "brms-port/results-merged.rds"))
k <- intersect(M$id, B$id)
then <- B$bucket[match(k, B$id)]; now <- M$bucket[match(k, M$id)]
kind <- M$kind[match(k, M$id)]
cat("other rows, 0.34.0 x 0.67.0:\n")
print(table(then[kind == "other"], now[kind == "other"]))
# ran at 0.34.0 in either pass, fails in both passes now
ran_then <- B$status_raw[match(k, B$id)] == "OK" |
  B$status_spell[match(k, B$id)] == "OK"
fails_now <- M$status_raw[match(k, M$id)] != "OK" &
  M$status_spell[match(k, M$id)] != "OK"
x <- k[ran_then & fails_now & kind %in% c("model", "post", "other")]
cat("\nran at 0.34.0 (raw or patched), fails raw and patched now:",
    length(x), "\n")
for (id in x) {
  cat(sprintf("%-26s %-5s %s -> %s | %s\n", id, M$kind[M$id == id],
              B$class[B$id == id], M$class[M$id == id],
              substr(M$msg_spell[M$id == id], 1, 100)))
}
