# Re-check: family constructors with arguments other than links, which
# family_call_of() would not carry into a stored bf() call.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
for (nm in getNamespaceExports("frmtmb")) {
  f <- tryCatch(get(nm, ns), error = function(e) NULL)
  if (!is.function(f)) next
  r <- tryCatch(suppressWarnings(suppressMessages(f())), error = function(e) NULL)
  if (!inherits(r, "frmtmb_family") && !inherits(r, "family")) next
  a <- setdiff(names(formals(f)), c("link", grep("^link_", names(formals(f)), value = TRUE)))
  if (length(a)) cat(nm, ":", a, "\n")
}
