# Compare two dumps of dev/spatial-compat-dump.R: every unordered pair's
# status and note. Rule 4 of a move: frm_compat() reports the same rows
# with the extension loaded as before.
#   Rscript dev/spatial-compat-compare.R <before.rds> <after.rds>
a <- commandArgs(trailingOnly = TRUE)
b <- readRDS(a[1]); f <- readRDS(a[2])
cat("before:", b$lib, b$version, paste(b$pkgs, collapse = ","), "\n")
cat("after: ", f$lib, f$version, paste(f$pkgs, collapse = ","), "\n")
cat("pairs", nrow(b$tab), "->", nrow(f$tab), "; features",
    nrow(b$features), "->", nrow(f$features), "\n")
fb <- b$features[order(b$features$name), c("name", "key", "kind")]
ff <- f$features[order(f$features$name), c("name", "key", "kind")]
rownames(fb) <- rownames(ff) <- NULL
cat("features identical (name, key, kind):", identical(fb, ff), "\n")
if (!identical(fb, ff)) {
  print(setdiff(fb$name, ff$name)); print(setdiff(ff$name, fb$name))
}
m <- merge(b$tab[, c("pair", "status", "note")],
           f$tab[, c("pair", "status", "note")], by = "pair", all = TRUE,
           suffixes = c(".b", ".a"))
only_b <- sum(is.na(m$status.a)); only_a <- sum(is.na(m$status.b))
both <- m[!is.na(m$status.a) & !is.na(m$status.b), ]
ds <- both[both$status.b != both$status.a, ]
dn <- both[both$status.b == both$status.a & both$note.b != both$note.a, ]
cat("pairs only before:", only_b, " only after:", only_a,
    " status differs:", nrow(ds), " note differs:", nrow(dn), "\n")
if (only_b) print(utils::head(m$pair[is.na(m$status.a)], 20))
if (only_a) print(utils::head(m$pair[is.na(m$status.b)], 20))
if (nrow(ds)) print(utils::head(ds[, c("pair", "status.b", "status.a")], 20))
if (nrow(dn)) print(utils::head(dn$pair, 20))
