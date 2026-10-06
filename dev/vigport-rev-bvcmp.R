# Reviewer: compare the reviewer's rerun of the hand translation on
# rellib-r5 with the lane's, row by row (vignette, label, ok).
#
#   Rscript dev/vigport-rev-bvcmp.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
a <- read.csv(file.path(root, "vigport-bv-out/r5/drift.csv"),
              stringsAsFactors = FALSE)
b <- read.csv(file.path(root, "vigport-rev-out/bv-r5/drift.csv"),
              stringsAsFactors = FALSE)
cat("lane rows", nrow(a), " reviewer rows", nrow(b), "\n")
key <- function(d) paste(d$vignette, d$label, ave(seq_len(nrow(d)),
                                                    d$vignette, d$label,
                                                    FUN = seq_along))
a$key <- key(a); b$key <- key(b)
cat("keys only lane:", length(setdiff(a$key, b$key)), " only reviewer:",
    length(setdiff(b$key, a$key)), "\n")
m <- merge(a[, c("key", "ok", "drift", "msg")], b[, c("key", "ok", "drift",
                                                       "msg")],
           by = "key", suffixes = c(".lane", ".rev"))
cat("ok at 0.67.0: lane", sum(a$ok), " reviewer", sum(b$ok), "\n")
d <- m[m$ok.lane != m$ok.rev, ]
cat("rows whose ok differs:", nrow(d), "\n")
for (i in seq_len(nrow(d))) cat(" ", d$key[i], d$ok.lane[i], d$ok.rev[i],
                                substr(d$msg.rev[i], 1, 120), "\n")
cat("\ndrift table, reviewer:\n")
b$path <- ifelse(grepl("^SAMPLE:", b$label), "SAMPLE", "ML")
print(table(b$drift, b$path))
cat("\nmessages differ on failing rows:",
    sum(!m$ok.rev & !m$ok.lane & m$msg.lane != m$msg.rev), "\n")
