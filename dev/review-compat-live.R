# Reviewer scratch: every sampled pair's live status, core-only and
# all-loaded, on the lane build; and the conditional live rows.
a <- readRDS("dev/review-compat-after.rds")
b <- readRDS("dev/review-compat-before.rds")
smp <- read.delim("dev/equiv-data/sample.tsv", stringsAsFactors = FALSE)
rec <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
look <- function(cp, x, y) {
  s <- cp$status[(cp$feature_a == x & cp$feature_b == y) |
                   (cp$feature_a == y & cp$feature_b == x)]
  if (length(s) != 1L) NA_character_ else s
}
smp$live_core <- mapply(function(x, y) look(a$core, x, y), smp$feature_a,
                        smp$feature_b)
smp$live_all <- mapply(function(x, y) look(a$all, x, y), smp$feature_a,
                       smp$feature_b)
smp$base_core <- mapply(function(x, y) look(b$core, x, y), smp$feature_a,
                        smp$feature_b)
smp$rec_verdict <- rec$verdict[match(smp$id, rec$id)]
smp$rec_status <- rec$status[match(smp$id, rec$id)]
cat("sampled:", nrow(smp), "\n")
cat("no core row (extension pairs):", sum(is.na(smp$live_core)), "\n")
print(smp[is.na(smp$live_core), c("id", "feature_a", "feature_b", "status",
                                  "live_all")])
cat("live_core != live_all:", sum(smp$live_core != smp$live_all,
                                  na.rm = TRUE), "\n")
cat("live (core) differs from stratum:",
    sum(smp$live_core != smp$status, na.rm = TRUE), "\n")
cat("live (lane) differs from base (0.80.0):",
    sum(smp$live_core != smp$base_core, na.rm = TRUE), "\n")
print(smp[which(smp$live_core != smp$base_core),
          c("id", "feature_a", "feature_b", "status", "base_core",
            "live_core")])
cat("\nlive != stratum:\n")
d <- smp[which(smp$live_core != smp$status), ]
print(d[, c("id", "feature_a", "feature_b", "status", "live_core",
            "rec_verdict")], row.names = FALSE)
cat("\nconditional live rows:\n")
cnd <- smp[which(smp$live_core == "conditional"), ]
print(cnd[, c("id", "feature_a", "feature_b", "status", "rec_verdict")],
      row.names = FALSE)
cat("\nstratum table:\n")
print(table(stratum = smp$status, live = smp$live_core, useNA = "ifany"))
