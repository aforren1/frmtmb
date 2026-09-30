# Punch round 2: which of the three families the round-1
# family_call_of() (its link-only comparison, reproduced here) wrote as
# a call.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
lk <- function(f) vapply(f[["links"]], function(l) {
  if (is.list(l)) l[["name"]] %||% NA_character_ else as.character(l)
}, "")
old <- function(fam) {
  nm <- fam[["family"]]
  ctor <- get(nm, envir = asNamespace("frmtmb"))
  links <- lk(fam)
  cl <- as.call(list(call("::", as.name("frmtmb"), as.name(nm))))
  for (dp in names(links)) {
    a <- if (dp == "mu") "link" else paste0("link_", dp)
    if (a %in% names(formals(ctor))) cl[[a]] <- links[[dp]]
  }
  back <- tryCatch(frmtmb:::as_frmtmb_family(eval(cl, baseenv())),
                   error = function(e) conditionMessage(e))
  if (is.character(back)) return(paste("eval error:", back))
  if (!identical(back[["family"]], nm) || !identical(lk(back), links)) {
    return("object")
  }
  deparse1(cl)
}
for (f in list(huber(k = 3), whittle(tapers = 4), cox(df = 6))) {
  cat(f$family, ":", old(frmtmb:::as_frmtmb_family(f)), "\n")
}
