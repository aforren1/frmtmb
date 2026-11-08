A <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-fx-band-after.rds")
B <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-fx-band-before.rds")
df <- do.call(rbind, lapply(names(A), function(id) data.frame(id = id, k = A[[id]]$k,
  basis = A[[id]]$basis, llA = A[[id]]$ll, llB = B[[id]]$ll, ref = A[[id]]$llref,
  wA = A[[id]]$warned, wB = B[[id]]$warned, errA = A[[id]]$err, errB = B[[id]]$err)))
cat("cells", nrow(df), " basis on (after):", sum(df$basis, na.rm = TRUE),
    " kappa 1e3..1e4:", sum(df$k >= 1e3 & df$k < 1e4), "\n")
cat("errors after", sum(df$errA), "before", sum(df$errB), "\n")
band <- df[df$basis %in% TRUE & df$k < 1e4, ]
cat("engaged at kappa < 1e4:", nrow(band), "\n")
d <- band$llA - band$llB
cat("|ll after - ll before| > 1e-6:", sum(abs(d) > 1e-6, na.rm = TRUE),
    " after higher:", sum(d > 1e-6, na.rm = TRUE), " after lower:", sum(d < -1e-6, na.rm = TRUE), "\n")
big <- band[abs(d) > 1e-6 & !is.na(d), ]
big$dAB <- big$llA - big$llB; big$dA_ref <- big$llA - big$ref; big$dB_ref <- big$llB - big$ref
print(big[order(-abs(big$dAB)), c("id", "k", "dAB", "dA_ref", "dB_ref", "wA", "wB")], row.names = FALSE, digits = 3)
cat("after vs reference, all cells: max short", min(df$llA - df$ref, na.rm = TRUE),
    " silent cells short > 1e-6:", sum(df$llA - df$ref < -1e-6 & !df$wA, na.rm = TRUE), "\n")
cat("before vs reference: max short", min(df$llB - df$ref, na.rm = TRUE),
    " silent cells short > 1e-6:", sum(df$llB - df$ref < -1e-6 & !df$wB, na.rm = TRUE), "\n")
cat("warned after", sum(df$wA), "before", sum(df$wB), "\n")
print(df[df$wA != df$wB, c("id", "k", "basis", "wA", "wB")], row.names = FALSE)
print(df[df$wA, c("id", "k", "basis", "llA", "ref", "wB")], row.names = FALSE, digits = 10)
