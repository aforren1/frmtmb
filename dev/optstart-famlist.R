# Lane optstart: every core family constructor, its dpars, links and
# init_dpars, to decide which families the least-squares start reaches.
#   Rscript dev/optstart-famlist.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
ex <- getNamespaceExports("frmtmb")
for (nm in sort(ex)) {
  f <- get(nm, ns)
  if (!is.function(f)) next
  fam <- tryCatch(suppressWarnings(suppressMessages(f())),
                  error = function(e) NULL)
  if (!inherits(fam, "frmtmb_family")) next
  lk <- vapply(fam[["links"]], function(l) l$name %||% "?", "")
  cat(sprintf("%-28s dpars %s | links %s | init %s | type %s\n", nm,
              paste(fam[["dpars"]], collapse = ","),
              paste(lk, collapse = ","),
              paste(names(fam[["init_dpars"]]), collapse = ","),
              fam[["type"]]))
}
