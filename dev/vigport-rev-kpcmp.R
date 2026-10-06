# Reviewer: keep-prior and need passes, lane against reviewer rerun.
#
#   Rscript dev/vigport-rev-kpcmp.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
rd <- function(d) {
  out <- list()
  for (f in list.files(d, "[.]rds$", full.names = TRUE)) {
    out <- c(out, readRDS(f))
  }
  out
}
for (m in c("results-keepprior", "results-need", "results", "results-spell")) {
  a <- rd(file.path(root, "vigport-port-out/r5", m))
  b <- rd(file.path(root, "vigport-rev-out/seeded", m))
  sa <- vapply(a, function(r) r$status, ""); sb <- vapply(b, function(r)
    r$status, "")
  ka <- vapply(a, function(r) r$kind, "")
  k <- intersect(names(sa), names(sb))
  cat(sprintf("%-18s lane %d rev %d common %d status differs %d\n", m,
              length(sa), length(sb), length(k), sum(sa[k] != sb[k])))
  for (kk in c("model", "post")) {
    sel <- k[ka[k] == kk & sb[k] != "SETUP-SKIP"]
    cat(sprintf("   %-5s rev OK %d of %d\n", kk, sum(sb[sel] == "OK"),
                length(sel)))
  }
}
