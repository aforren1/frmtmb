# frm_compat() as a session reports it, keyed by the unordered pair, so
# that the base session (core plus frmtmb.spline, rellib-r11) and the
# lane's session (core plus frmtmb.survival) can be compared row for row.
#
#   Rscript dev/splitsurv-compat-dump.R <before|after> <out.rds>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
out <- a[2]
LIB <- "C:/Users/adf44/source/r/wt-splitsurv-lib"
REL <- "C:/Users/adf44/source/r/rellib-r11"
USR <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "after")) LIB, REL, USR))
pkg <- if (identical(arm, "after")) "frmtmb.survival" else "frmtmb.spline"
suppressMessages(library(pkg, character.only = TRUE))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
tab <- frm_compat()
key <- ifelse(tab$feature_a < tab$feature_b,
              paste(tab$feature_a, tab$feature_b, sep = " x "),
              paste(tab$feature_b, tab$feature_a, sep = " x "))
tab$key <- key
ft <- frm_compat_features()
watch <- c("cox", "royston_parmar", "rp_floored", "frm_curve",
           "frm_curve_contrast")
sel <- tab$feature_a %in% watch | tab$feature_b %in% watch
cat("pairs:", nrow(tab), " features:", nrow(ft), " watched pairs:",
    sum(sel), "\n")
for (w in watch) {
  cat(sprintf("%-20s kind=%-8s pairs=%d\n", w,
              paste(ft$kind[ft$name == w], collapse = ","),
              sum(tab$feature_a == w | tab$feature_b == w)))
}
saveRDS(list(tab = tab, features = ft, unresolved = attr(tab, "unresolved")),
        out)
