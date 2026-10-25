## Base vs lane prior tables of dev/review-brmsnames-fr1.R and -fr1vp.R:
## the same rows (base coef in brms's spelling), only the order changed.
.libPaths(c("C:/Users/adf44/source/r/wt-brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pcb <- getFromNamespace("prior_coef_brms", "frmtmb")
cols <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar",
          "lb", "ub", "source")
canon <- function(t, rename) {
  t <- as.data.frame(t)[, intersect(cols, names(as.data.frame(t))),
                        drop = FALSE]
  # the base table's coef in brms's spelling, its alias dropped
  if (rename) {
    cf <- gsub(".frm_", "", t$coef, fixed = TRUE)
    t$coef <- ifelse(nzchar(cf), pcb(cf), cf)
  }
  t[] <- lapply(t, as.character)
  r <- do.call(paste, c(t, sep = "|"))
  list(sorted = sort(r), raw = r)
}
for (stem in c("fr1", "fr1vp")) {
  b <- readRDS(sprintf("dev/review-brmsnames-out/%s-tables-base.rds", stem))
  l <- readRDS(sprintf("dev/review-brmsnames-out/%s-tables-lane.rds", stem))
  for (nm in names(b)) for (w in c("dp", "vp")) {
    if (is.character(b[[nm]][[w]]) || is.character(l[[nm]][[w]])) next
    cb <- canon(b[[nm]][[w]], TRUE)
    cl <- canon(l[[nm]][[w]], FALSE)
    cat(sprintf("%-6s %-2s %-66s rows %d/%d same rows: %s; same order: %s\n",
                stem, w, substr(nm, 1, 66), length(cb$raw), length(cl$raw),
                identical(cb$sorted, cl$sorted),
                identical(cb$raw, cl$raw)))
  }
}
