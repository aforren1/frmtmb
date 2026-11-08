# Lane bigx: the full warning of the grid's cells at 1e8 that warn
# (dev/bigx-out/grid-after5.tsv), to read the clause after "A Newton".
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
cells <- list(c("poisson", TRUE, 1), c("gaussian", FALSE, 3),
              c("Gamma", TRUE, 3))
for (cl in cells) {
  re <- as.logical(cl[2])
  d <- bigx_data(cl[1], 1e8, as.integer(cl[3]), re = re)
  rhs <- if (re) " + (1 | g)" else ""
  w <- character(0)
  withCallingHandlers(
    frm(bf(as.formula(paste0("y ~ w", rhs))), family = bigx_family(cl[1]),
        data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  cat(cl, ":\n")
  writeLines(strwrap(w, 76, prefix = "  "))
}
