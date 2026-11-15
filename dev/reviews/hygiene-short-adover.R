# Short check of lane hygiene's fix round, item 5: where RTMB attaches
# and detaches "AD-overloads", and which library hmmTMB comes from.
#   Rscript dev/reviews/hygiene-short-adover.R
.libPaths(c("C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("RTMB")
for (n in ls(ns, all.names = TRUE)) {
  f <- get(n, ns)
  if (!is.function(f)) next
  d <- deparse(f)
  if (!any(grepl("AD-overloads", d, fixed = TRUE))) next
  cat("== RTMB:::", n, "\n", sep = "")
  i <- grep("AD-overloads|attach|detach|on.exit", d)
  cat(d[sort(unique(i))], sep = "\n")
}
cat("hmmTMB from", find.package("hmmTMB"), "\n")
cat("RTMB", as.character(packageVersion("RTMB")), "from",
    find.package("RTMB"), "\n")
