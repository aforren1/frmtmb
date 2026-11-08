# The sampled equivalence pairs whose spec carries a prior, with their
# recorded verdicts: the pairs the split can reach (lane brief, "Rerun
# the sampled pairs that have a prior").
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R")) sys.source(f, envir = env)
smp <- read.delim("dev/equiv-data/sample.tsv")
res <- read.delim("dev/equiv-data/results.tsv")
keys <- paste(smp$feature_a, smp$feature_b, sep = " x ")
has <- vapply(keys, function(k) {
  fn <- env$eq_specs[[k]]
  if (is.null(fn)) return(NA)
  sp <- tryCatch(fn(), error = function(e) NULL)
  !is.null(sp) && !is.null(sp$prior)
}, NA)
out <- data.frame(id = smp$id, key = keys, prior = has)
out <- merge(out, res[, c("id", "verdict", "level", "l1_const", "l4_z")],
             by = "id", all.x = TRUE)
out <- out[which(out$prior), ]
cat("pairs with a prior:", nrow(out), "of", nrow(smp), "\n")
write.table(out, "dev/loglik-r80-equiv-ids.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
print(out[, c("id", "key", "verdict", "level")], row.names = FALSE)
