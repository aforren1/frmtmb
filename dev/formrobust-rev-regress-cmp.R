# Reviewer: compare regress-base.rds and regress-lane.rds item by item
# with identical(); print every difference and every error.
.libPaths(c("C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-log/"
b <- readRDS(paste0(p, "regress-base.rds"))
l <- readRDS(paste0(p, "regress-lane.rds"))
items <- c("par", "fn", "gr", "fn2", "gr2", "estimates", "logLik", "fixef",
           "fitted", "fitted_nd", "predict", "predict_nd", "residuals",
           "simulate", "ce", "emm", "emm_x", "emm_ep", "default_prior",
           "data_names")
nident <- 0L; ndiff <- 0L; nerr_both <- 0L
for (m in names(b)) {
  for (it in items) {
    x <- b[[m]][[it]]; y <- l[[m]][[it]]
    ex <- inherits(x, "rev_err"); ey <- inherits(y, "rev_err")
    if (ex && ey && identical(x, y)) {
      nerr_both <- nerr_both + 1L
      cat(sprintf("%-14s %-13s both error: %s\n", m, it,
                  substr(x, 1, 110)))
      next
    }
    if (identical(x, y)) { nident <- nident + 1L; next }
    ndiff <- ndiff + 1L
    cat(sprintf("%-14s %-13s DIFFERS", m, it))
    if (ex || ey) {
      cat(sprintf("\n    base: %s\n    lane: %s\n",
                  substr(if (ex) x else "<value>", 1, 200),
                  substr(if (ey) y else "<value>", 1, 200)))
    } else {
      ae <- all.equal(x, y)
      cat(": ", paste(head(ae, 3), collapse = " | "), "\n")
    }
  }
}
cat(sprintf("\nmodels %d, items %d: identical %d, both-error %d, differ %d\n",
            length(b), length(b) * length(items), nident, nerr_both, ndiff))
