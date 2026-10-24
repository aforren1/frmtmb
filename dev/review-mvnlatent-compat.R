# Reviewer's independent compat comparison.
#   Rscript dev/review-mvnlatent-compat.R <lib or base> <out.rds> [pkg ...]
#   Rscript dev/review-mvnlatent-compat.R compare <a.rds> <b.rds> [drop]
a <- commandArgs(trailingOnly = TRUE)
if (a[1] != "compare") {
  base <- c("C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6")
  .libPaths(if (identical(a[1], "base")) base else c(a[1], base))
  suppressMessages(library(frmtmb))
  for (p in a[-(1:2)]) suppressMessages(loadNamespace(p))
  cat("frmtmb", format(packageVersion("frmtmb")),
      dirname(find.package("frmtmb")), "\n")
  for (p in a[-(1:2)]) cat(p, dirname(find.package(p)), "\n")
  saveRDS(list(tb = as.data.frame(frm_compat()),
               ft = as.data.frame(frm_compat_features())), a[2])
  quit(save = "no")
}
A <- readRDS(a[2])
B <- readRDS(a[3])
cat("columns tb equal:", identical(names(A$tb), names(B$tb)), "\n")
cat("columns:", paste(names(A$tb), collapse = ","), "\n")
key <- function(t) {
  fa <- pmin(t$feature_a, t$feature_b)
  fb <- pmax(t$feature_a, t$feature_b)
  ka <- ifelse(t$feature_a <= t$feature_b, t$kind_a, t$kind_b)
  kb <- ifelse(t$feature_a <= t$feature_b, t$kind_b, t$kind_a)
  rest <- t[setdiff(names(t), c("feature_a", "feature_b", "kind_a",
                                "kind_b"))]
  do.call(paste, c(list(fa, ka, fb, kb), rest, sep = "\r"))
}
ta <- A$tb
if (length(a) > 3 && a[4] == "drop") {
  ta <- ta[ta$feature_a != "mixture_mvn" & ta$feature_b != "mixture_mvn", ]
  A$ft <- A$ft[A$ft$name != "mixture_mvn", ]
}
ka <- key(ta)
kb <- key(B$tb)
cat("rows", length(ka), length(kb), " dup keys", anyDuplicated(ka),
    anyDuplicated(kb), "\n")
cat("only A", length(setdiff(ka, kb)), " only B", length(setdiff(kb, ka)),
    "\n")
print(utils::head(setdiff(ka, kb)))
print(utils::head(setdiff(kb, ka)))
fk <- function(f) do.call(paste, c(f, sep = "\r"))
cat("vocab", nrow(A$ft), nrow(B$ft), " only A",
    length(setdiff(fk(A$ft), fk(B$ft))), " only B",
    length(setdiff(fk(B$ft), fk(A$ft))), "\n")
print(utils::head(setdiff(fk(A$ft), fk(B$ft))))
print(utils::head(setdiff(fk(B$ft), fk(A$ft))))
