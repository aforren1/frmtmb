# Reviewer: does the stratum status of a core pair change when the
# extensions load? The dump takes status from the all-loaded table.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
t0 <- as.data.frame(frm_compat())
cat("core-only:", nrow(t0), "\n")
print(table(t0$status))
cp <- read.delim("dev/equiv-frame-out/compat.tsv")
core <- cp[cp$in_core_only, ]
cat("in_core_only rows:", nrow(core), "\n")
print(table(core$status))
k0 <- paste(t0$feature_a, t0$feature_b)
k1 <- paste(core$feature_a, core$feature_b)
k1r <- paste(core$feature_b, core$feature_a)
m <- match(k1, k0)
m[is.na(m)] <- match(k1r[is.na(m)], k0)
cat("unmatched:", sum(is.na(m)), "\n")
d <- core$status != t0$status[m]
cat("status differs core-only vs all-loaded:", sum(d), "\n")
if (any(d)) print(head(cbind(core[d, c("feature_a", "feature_b", "status")],
                             core_status = t0$status[m][d]), 40))
smp <- read.delim("dev/equiv-data/sample.tsv")
ks <- paste(smp$feature_a, smp$feature_b)
ms <- match(ks, k0)
ms[is.na(ms)] <- match(paste(smp$feature_b, smp$feature_a)[is.na(ms)], k0)
ds <- smp$status != t0$status[ms]
cat("sampled pairs whose core-only status differs:", sum(ds), "\n")
if (any(ds)) print(cbind(smp[ds, c("id", "feature_a", "feature_b", "status")],
                         core_status = t0$status[ms][ds]))
# duplicate pairs in either orientation?
cat("dup unordered pairs in core:",
    sum(duplicated(apply(core[, c("feature_a", "feature_b")], 1,
                         function(x) paste(sort(x), collapse = " ")))), "\n")
