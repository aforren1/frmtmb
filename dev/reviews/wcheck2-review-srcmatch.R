# Review of lane wcheck2: the installed lane library holds the tree's code
.libPaths(c("C:/Users/adf44/source/r/wt-wcheck2-lib",
            "C:/Users/adf44/source/r/rellib-r21", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
e <- new.env()
for (f in c("R/weights-unbounded.R")) sys.source(f, envir = e, keep.source = FALSE)
src <- new.env()
exprs <- parse("R/fit.R", keep.source = FALSE)
for (x in exprs) if (is.call(x) && identical(x[[1]], as.name("<-")) &&
                     as.character(x[[2]]) %in% c("check_re_structure",
                                                  "re_check_term")) eval(x, src)
for (nm in c(ls(e), ls(src))) {
  a <- get(nm, ns); b <- if (exists(nm, e, inherits = FALSE)) get(nm, e) else get(nm, src)
  same <- identical(deparse(a, control = NULL), deparse(b, control = NULL))
  if (!same) cat("DIFFERS:", nm, "\n")
}
cat("compared", length(ls(e)) + length(ls(src)), "objects\n")
