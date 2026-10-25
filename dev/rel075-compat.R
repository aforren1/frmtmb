# The merged compatibility table against 0.74.0, with every extension
# loaded (rule 4 of a move, dev/core-boundary.md). The base is
# rellib-r12 (frmtmb 0.74.0 and its eight extensions), the merged build
# rellib-r13 (frmtmb 0.75.0 and its eight). Adapted from
# dev/rel074-compat.R. This round moved nothing, but three lanes change
# the table on purpose:
#   bhaz      bhaz() is frmtmb.survival's option term of cox(); its
#             rules, and the notes of core's four generic rules for an
#             addition term (kind:aterm), which now say what a data term
#             and an option term each do;
#   sarfcor   fcor() is core's and sar() frmtmb.spatial's, residual
#             correlation terms with their own rules;
#   smallfix, brmsnames: no compat rule (their records say so).
# So every pair whose status or note differs must have bhaz(), fcor()
# or sar() on one side, or (a note only) an addition term on one side.
# Anything else fails the verdict.
#   Rscript dev/rel075-compat.R      (after the two dumps:)
#   Rscript dev/rel075-compat-dump.R dev/rel075-log/compat-base.rds \
#     C:/Users/adf44/source/r/rellib-r12 <the eight>
#   Rscript dev/rel075-compat-dump.R dev/rel075-log/compat-merged.rds \
#     C:/Users/adf44/source/r/rellib-r13 <the eight>
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel075")
tab_of <- function(x) {
  t <- x$tab
  t$pair <- ifelse(t$feature_a < t$feature_b,
                   paste(t$feature_a, t$feature_b, sep = " x "),
                   paste(t$feature_b, t$feature_a, sep = " x "))
  t$orient <- paste(t$feature_a, t$feature_b, sep = " x ")
  stopifnot(!anyDuplicated(t$pair))
  t[, c("pair", "orient", "feature_a", "feature_b", "status", "note")]
}
base <- readRDS("dev/rel075-log/compat-base.rds")
mrg <- readRDS("dev/rel075-log/compat-merged.rds")
cat("base:  ", paste(base$libs, collapse = "\n        "), "\n")
cat("merged:", paste(mrg$libs, collapse = "\n        "), "\n")
B <- tab_of(base); M <- tab_of(mrg)
m <- merge(B, M, by = "pair", all = TRUE, suffixes = c(".b", ".a"))
only_b <- m$pair[is.na(m$status.a)]
only_a <- m$pair[is.na(m$status.b)]
both <- m[!is.na(m$status.a) & !is.na(m$status.b), ]
cat("pairs", nrow(B), "->", nrow(M), "| only base", length(only_b),
    "| only merged", length(only_a), "\n")
fb <- base$features[, c("name", "kind")]
fm <- mrg$features[, c("name", "kind")]
fx <- merge(fb, fm, by = "name", all = TRUE, suffixes = c(".b", ".a"))
fd <- fx[is.na(fx$kind.b) | is.na(fx$kind.a) | fx$kind.b != fx$kind.a, ]
cat("features", nrow(fb), "->", nrow(fm), "; added, removed or of",
    "another kind:\n")
print(fd, row.names = FALSE)
new_feat <- c("bhaz()", "fcor()", "sar()")
aterms <- fm$name[fm$kind == "aterm"]
has <- function(p, set) {
  vapply(strsplit(p, " x ", fixed = TRUE), function(v) any(v %in% set), NA)
}
st <- both[both$status.b != both$status.a, ]
nt <- both[both$note.b != both$note.a & both$status.b == both$status.a, ]
cat("\nstatus differs:", nrow(st), "| with bhaz(), fcor() or sar():",
    sum(has(st$pair, new_feat)), "\n")
if (nrow(st)) print(table(paste(st$status.b, "->", st$status.a),
            ifelse(has(st$pair, "bhaz()"), "bhaz()",
                   ifelse(has(st$pair, "fcor()"), "fcor()",
                          ifelse(has(st$pair, "sar()"), "sar()", "OTHER")))))
cat("\nnote differs (status the same):", nrow(nt),
    "| with bhaz(), fcor() or sar():", sum(has(nt$pair, new_feat)),
    "| other, with an addition term:",
    sum(!has(nt$pair, new_feat) & has(nt$pair, aterms)), "\n")
other_note <- nt[!has(nt$pair, new_feat) & !has(nt$pair, aterms), ]
cat("note differs on a pair with neither:", nrow(other_note), "\n")
for (i in seq_len(min(nrow(other_note), 20))) {
  cat("  ", other_note$pair[i], "\n    base:  ",
      substr(other_note$note.b[i], 1, 160), "\n    merged:",
      substr(other_note$note.a[i], 1, 160), "\n")
}
# the generic aterm notes: the rewording bhaz made, on every pair it
# reaches (one example)
ga <- nt[!has(nt$pair, new_feat) & has(nt$pair, aterms), ]
if (nrow(ga)) {
  cat("example of an addition-term note:", ga$pair[1], "\n    base:  ",
      substr(ga$note.b[1], 1, 200), "\n    merged:",
      substr(ga$note.a[1], 1, 200), "\n")
  cat("distinct (base note, merged note) pairs among them:",
      nrow(unique(ga[, c("note.b", "note.a")])), "\n")
}
cat("pairs only in the merged table, by new feature:\n")
print(table(ifelse(has(only_a, "bhaz()"), "bhaz()",
                   ifelse(has(only_a, "fcor()"), "fcor()",
                          ifelse(has(only_a, "sar()"), "sar()", "OTHER")))))
fl <- both$pair[both$orient.b != both$orient.a]
cat("orientation flips (no status or note meaning):", length(fl), "\n")
ok <- length(only_b) == 0 && all(has(only_a, new_feat)) &&
  all(fd$name %in% new_feat) && all(has(st$pair, new_feat)) &&
  nrow(other_note) == 0
cat("COMPAT VERDICT:", if (ok) "PASS" else "FAIL", "\n")
