# Review lane osa, short check: the message frm() gives for cens() and
# trunc() on every core family, so the two arms can be diffed. Data are
# chosen to be valid for most families (0/1 counts, values in (0, 1));
# what matters is only whether the message changes between arms.
#   Rscript dev/review-osa-aterms.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-aterms", arm, "\n")
set.seed(1)
d <- data.frame(x = rnorm(40), cc = 0, k = 3L,
                yd = rep(0:1, 20), yc = runif(40, 0.1, 0.9))
reg <- frmtmb:::family_registry
for (nm in sort(names(reg))) {
  fam <- tryCatch(reg[[nm]](), error = function(e) NULL)
  if (is.null(fam) || !is.list(fam)) next
  y <- if (identical(fam$type, "discrete") || identical(fam$type, "ordinal")) "yd" else "yc"
  for (at in c("cens(cc)", "trunc(lb = 0)")) {
    f <- stats::as.formula(paste0(y, " | ", at, " ~ x"))
    m <- tryCatch({
      suppressWarnings(frm(f, family = fam, data = d))
      "fits"
    }, error = function(e) gsub("\n", " ", conditionMessage(e)))
    cat(sprintf("%-28s %-14s %s\n", nm, at, substr(m, 1, 150)))
  }
}
cat("END review-osa-aterms\n")
