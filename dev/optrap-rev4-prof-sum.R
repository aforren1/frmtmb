# Reviewer of lane optrap, consolidation check: newton_polish() calls of
# the core suite on build 32 (dev/optrap-rev4-profile.R logs).
fs <- Sys.glob("dev/optrap-rev4-log/prof-core/*.txt")
X <- do.call(rbind, lapply(fs, function(f) {
  read.delim(f, header = FALSE, quote = "", stringsAsFactors = FALSE,
             col.names = c("k", "tag", "re", "ill", "np", "code", "steps",
                           "before", "after", "secs"))
}))
X$d <- X$after - X$before
cat("calls", nrow(X), "in", length(unique(X$tag)), "test files\n")
print(table(re = X$re, ill = X$ill, stepped = X$steps > 0))
r <- X[X$re == 1 & X$ill, ]
cat(sprintf("random effects, ill: %d calls in %d files | stepped %d | objective change min %.3g max %.3g | secs max %.3f, total %.2f\n",
            nrow(r), length(unique(r$tag)), sum(r$steps > 0), min(r$d),
            max(r$d), max(r$secs), sum(r$secs)))
print(table(r$tag, r$steps > 0))
f <- X[X$re == 0, ]
cat(sprintf("no random effects: %d calls | stepped %d | change min %.3g max %.3g | secs max %.3f\n",
            nrow(f), sum(f$steps > 0), min(f$d), max(f$d), max(f$secs)))
cat("any call with a higher objective after:", sum(X$d > 0), "\n")
ok <- is.finite(X$d)
cat("finite objective changes:", sum(ok), "of", nrow(X), "| higher after:",
    sum(X$d[ok] > 0), "| stepped calls' change range:",
    paste(signif(range(X$d[ok & X$steps > 0]), 3), collapse = " to "), "\n")
print(X[!ok, c("tag", "re", "steps", "before", "after")][1:5, ])
