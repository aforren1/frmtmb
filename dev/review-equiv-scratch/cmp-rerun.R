# Reviewer: compare the warm rerun (dev/review-equiv-out) with the
# lane's reported results (dev/equiv-out), column by column.
rd <- function(p) read.delim(p, stringsAsFactors = FALSE)
ids <- sub("[.]tsv$", "", list.files("dev/review-equiv-out2", "[.]tsv$"))
cols <- c("verdict", "level", "agree", "frm_ok", "brms_ok", "l1_ok",
          "l2_ok", "l4_ok", "proposal")
num <- c("l1_const", "l1_grad", "l2_metric")
out <- list()
for (id in ids) {
  a <- rd(file.path("dev/equiv-out", paste0(id, ".tsv")))
  b <- rd(file.path("dev/review-equiv-out2", paste0(id, ".tsv")))
  same <- vapply(cols, function(k) identical(as.character(a[[k]]),
                                             as.character(b[[k]])), NA)
  dn <- vapply(num, function(k) {
    x <- as.numeric(a[[k]]); y <- as.numeric(b[[k]])
    if (is.na(x) && is.na(y)) 0 else abs(x - y)
  }, 0)
  out[[id]] <- data.frame(id = id, verdict = b$verdict, level = b$level,
                          same_cols = all(same),
                          differ = paste(cols[!same], collapse = ","),
                          d_const = signif(dn[1], 2), d_grad = signif(dn[2], 2),
                          d_l2 = signif(dn[3], 2),
                          l4z_lane = a$l4_z, l4z_rerun = b$l4_z,
                          sec_lane = a$seconds, sec_rerun = b$seconds)
}
res <- do.call(rbind, out)
options(width = 220)
print(res, row.names = FALSE)
cat("pairs:", nrow(res), " identical verdict columns:", sum(res$same_cols),
    "\n")
cat("rerun seconds sum:", sum(res$sec_rerun), " lane seconds sum:",
    sum(res$sec_lane), "\n")
