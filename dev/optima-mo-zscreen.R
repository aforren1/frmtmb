# Lane optima, item 1: could mo_search() skip a term whose coefficient
# is far from zero? For each seed of dev/optima-mo-study.R (and the
# plain ls ~ mo(income)), the chart-only fit (mo_search() disabled) and
# each mo() coefficient's z: the estimate over its standard error from
# fixef(), and over the conditional one, 1 / sqrt(H_bb) from the exact
# Hessian. Merged afterwards with the gains of the full search.
#   Rscript dev/optima-mo-zscreen.R int|main seeds out.tsv
args <- commandArgs(trailingOnly = TRUE)
form <- args[1]
seeds <- eval(parse(text = args[2]))
out <- args[3]
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
real <- utils::getFromNamespace("mo_search", "frmtmb")
src <- readLines("dev/optima-mo-study.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^# exact maximum", src) - 1L)]))
fo <- if (form == "int") ls ~ mo(income) * age else ls ~ mo(income)
rows <- list()
for (s in seeds) {
  d <- mk(s)
  utils::assignInNamespace("mo_search", function(obj, opt, ...) opt,
                           "frmtmb")
  f <- suppressWarnings(frm(fo, data = d))
  utils::assignInNamespace("mo_search", real, "frmtmb")
  g <- suppressWarnings(frm(fo, data = d))
  fe <- fixef(f)
  H <- f$obj$he(f$opt$par)
  terms <- frmtmb:::mo_search_terms(f$obj, f$frame)
  for (j in seq_along(terms)) {
    tm <- terms[[j]]
    b <- f$opt$par[tm$b]
    w <- frmtmb:::mo_simplex(f$opt$par[tm$z])
    rows[[length(rows) + 1L]] <- data.frame(
      seed = s, term = j, b = b,
      z_marg = b / fe[grep("^moincome", rownames(fe))[j], 2],
      z_cond = b * sqrt(H[tm$b, tm$b]), min_w = min(w),
      ll_chart = as.numeric(logLik(f)), ll_full = as.numeric(logLik(g)))
  }
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("rows", nrow(X), "\n")
