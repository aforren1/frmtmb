# Short check: the 53 proposals have no duplicate pair (either order),
# their status_now matches the pair table, and none contradicts a
# sampled verdict or the nl sweep.
p <- read.delim("dev/equiv-data/compat-proposals.tsv", stringsAsFactors = FALSE)
r <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
s <- read.delim("dev/equiv-data/nl-sweep.tsv", stringsAsFactors = FALSE)
cp <- read.delim("dev/equiv-frame-out/compat.tsv")
cp <- cp[cp$in_core_only, ]
key <- function(a, b) apply(cbind(a, b), 1, function(x)
  paste(sort(x), collapse = " x "))
p$k <- key(p$feature_a, p$feature_b)
cat("proposals:", nrow(p), " unique pairs:", length(unique(p$k)), "\n")
cp$k <- key(cp$feature_a, cp$feature_b)
st <- cp$status[match(p$k, cp$k)]
cat("status_now matches table:", sum(st == p$status_now), "of", nrow(p),
    "\n")
r$k <- key(r$feature_a, r$feature_b)
both <- intersect(p$k, r$k)
cat("proposals on sampled pairs:", length(both), "\n")
# sweep rows also present as sampled pairs
sk <- key(rep("nl", nrow(s)), s$family)
cat("sampled nl x family pairs:", paste(intersect(sk, r$k), collapse = "; "),
    "\n")
x <- merge(p[p$k %in% r$k, c("k", "proposal", "verdict")],
           r[, c("k", "verdict", "agree")], by = "k",
           suffixes = c(".prop", ".res"))
print(x[x$verdict.prop != x$verdict.res, ])
cat("sweep proposal vs proposals table mismatches:\n")
m <- match(key(rep("nl", nrow(s)), s$family), p$k)
print(data.frame(family = s$family, sweep = s$proposal,
                 table = p$proposal[m])[!is.na(m) &
                   s$proposal != p$proposal[m], ])
cat("honored but no proposal:", paste(s$family[s$outcome == "nl honored" &
  is.na(m) & s$status_now == "refused"], collapse = ", "), "\n")
